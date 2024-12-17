#pragma once

#include "general_utility.h"

// TODO: clean.up in another header
std::array<std::array<float, 2>, 2> get_intercept(TF1 *pol1, TF1 *pol2) { return utility::get_intercept(pol1->GetParameter(1), pol1->GetParError(1), pol1->GetParameter(0), pol1->GetParError(0), pol2->GetParameter(1), pol2->GetParError(1), pol2->GetParameter(0), pol2->GetParError(0)); }
std::array<std::array<float, 2>, 2> get_intercept_x(TF1 *pol1) { return utility::get_intercept_x(pol1->GetParameter(1), pol1->GetParError(1), pol1->GetParameter(0), pol1->GetParError(0)); }

namespace breakdown_voltage
{
    //  Basic method to find a first Vbd guess
    std::array<double, 2> find_vbd_guess(TGraphErrors *gLogTarget, float threshold = 0.35);

    //  Implemented methods to measure Vbd
    std::array<float, 2> measure_breakdown_0(TGraphErrors *gTarget, std::string image_folder = "", float threshold_guess = 0.75, bool erf_mediate = true, float min_before_fit = 5., float max_before_fit = 0.2, float min_after_fit = 0.2, float max_after_fit = 1.);
    std::array<float, 2> measure_breakdown_1(TGraphErrors *gTarget, std::string image_folder = "", float threshold_guess = 0.75);
    std::array<float, 2> measure_breakdown_2(TGraphErrors *gTarget, std::string image_folder = "", float threshold_guess = 0.75);
    std::array<float, 2> measure_breakdown_3(TGraphErrors *gTarget, std::string image_folder = "", float threshold_guess = 0.75);
    std::array<float, 2> measure_breakdown_4(TGraphErrors *gTarget, std::string image_folder = "", float threshold_guess = 0.75);
    std::array<float, 2> measure_breakdown_5(TGraphErrors *gTarget, std::string image_folder = "", float threshold_guess = 0.75);
}

std::array<double, 2> breakdown_voltage::find_vbd_guess(TGraphErrors *gLogTarget, float threshold)
{
    auto vbd_guess = 0.;
    auto mean_baseline = gLogTarget->GetY()[0];
    auto mean_contributors = 1;
    for (int iPnt = 1; iPnt < gLogTarget->GetN(); iPnt++)
    {
        auto current_Y = gLogTarget->GetY()[iPnt];
        if (fabs(mean_baseline - current_Y) > threshold)
        {
            vbd_guess = gLogTarget->GetX()[iPnt - 1];
            break;
        }
        mean_baseline *= mean_contributors;
        mean_baseline += current_Y;
        mean_contributors++;
        mean_baseline /= mean_contributors;
    }
    return {vbd_guess, mean_baseline};
}

std::array<float, 2> breakdown_voltage::measure_breakdown_0(TGraphErrors *gTarget, std::string image_folder = "", float threshold_guess = 0.75, bool erf_mediate = true, float min_before_fit = 5., float max_before_fit = 0.2, float min_after_fit = 0.2, float max_after_fit = 1.)
{
    //  Method 0:
    //  Linear fit before and after guess of breakdown voltage, take the intersection of the two as the Vbd value

    //  Create result container
    std::array<float, 2> result = {-1, -1};

    //  Make the log graph
    auto log_graph = graphutils::log(gTarget);

    //  Calculate a guess for the Vbd and have the measured baseline value
    auto guess_result = breakdown_voltage::find_vbd_guess(log_graph, threshold_guess);
    auto vbd_guess = guess_result[0];
    auto mean_baseline = guess_result[1];

    //  Define the functions to find the Vbd
    TF1 *pol1_before = new TF1("pol1_before", "[0]+x*[1]");
    TF1 *pol1_after = new TF1("pol1_after", "[0]+x*[1]");

    //  Fit before and after the Vbd guess
    //  --- before fit
    pol1_before->SetParameter(0, mean_baseline);
    pol1_before->SetParameter(1, 0.);
    log_graph->Fit(pol1_before, "RQ", "", vbd_guess - min_before_fit, vbd_guess - max_before_fit);
    //  --- after fit
    pol1_after->SetParameter(0, -100);
    pol1_after->SetParameter(1, 4);
    log_graph->Fit(pol1_after, "RQ", "", vbd_guess + min_after_fit, vbd_guess + max_after_fit);

    //  Intercept of the two lines
    auto intercept = get_intercept(pol1_before, pol1_after);
    result[0] = intercept[0][0];
    result[1] = intercept[0][1];

    TF1 *ffull;
    //  Mediate with erf if requested
    if (erf_mediate)
    {
        ffull = new TF1("ffull", "((0.5-0.5*TMath::Erf(((x-[0])/[1])))*([2]+x*[3])+(0.5+0.5*TMath::Erf(((x-[0])/[1])))*([4]+x*[5]))");
        ffull->SetParameter(0, intercept[0][0]);
        ffull->SetParLimits(0, intercept[0][0] - max_before_fit, intercept[0][0] + min_after_fit);
        ffull->SetParameter(1, 0.0000015);
        ffull->SetParLimits(1, 0.000001, 0.01);
        ffull->FixParameter(2, pol1_before->GetParameter(0));
        ffull->FixParameter(3, pol1_before->GetParameter(1));
        ffull->FixParameter(4, pol1_after->GetParameter(0));
        ffull->FixParameter(5, pol1_after->GetParameter(1));
        log_graph->Fit(ffull, "RQ", "SAME", result[0] - 4, result[0] + 3.);
        ffull->SetParLimits(0, intercept[0][0] + ffull->GetParameter(1), intercept[0][0] + 0.2);
        log_graph->Fit(ffull, "RQ", "SAME", result[0] - 4, result[0] + +3.);
        result[0] = ffull->GetParameter(0) - ffull->GetParameter(1);
        result[1] = sqrt(ffull->GetParError(0) * ffull->GetParError(0) + ffull->GetParameter(1) * ffull->GetParameter(1));
    }

    if (!image_folder.empty())
    {
        //  Create Canvas with result
        TCanvas *c1 = new TCanvas("c1", "c1", 600, 500);
        gPad->SetMargin(0.15, 0.05, 0.10, 0.05);

        //  Set-up graph graphics
        log_graph->GetXaxis()->SetTitle("bias voltage (V)");
        log_graph->GetYaxis()->SetTitle("log of current");
        log_graph->SetMarkerStyle(20);
        log_graph->SetMarkerColor(kAzure - 3);
        log_graph->GetXaxis()->SetRangeUser(result[0] - 6, result[0] + 3);
        log_graph->GetYaxis()->SetRangeUser(pol1_before->Eval(result[0] - 6) - 1, pol1_after->Eval(result[0] + 3) + 1);
        log_graph->Draw("ALPE");

        //  Plot all fits
        pol1_before->SetRange(result[0] - 5, result[0] + 5);
        pol1_before->SetLineColor(kBlue);
        pol1_before->DrawCopy("SAME");
        pol1_after->SetRange(result[0] - 5, result[0] + 5);
        pol1_after->SetLineColor(kGreen);
        pol1_after->DrawCopy("SAME");
        ffull->SetRange(result[0] - 5, result[0] + 5);
        ffull->SetLineColor(kRed);
        ffull->DrawCopy("SAME");

        //  Write on canvas results
        TLatex *l1 = new TLatex();
        l1->DrawLatexNDC(0.45, 0.90, Form("V_{bd} (V)"));
        l1->DrawLatexNDC(0.18, 0.85, Form("guess"));
        l1->DrawLatexNDC(0.35, 0.85, Form("= %.2f ", vbd_guess));
        l1->DrawLatexNDC(0.18, 0.80, Form("intercept"));
        l1->DrawLatexNDC(0.35, 0.80, Form("= %.2f #pm %.2f", intercept[0][0], intercept[0][1]));
        if (erf_mediate)
        {
            l1->DrawLatexNDC(0.18, 0.75, Form("fit"));
            l1->DrawLatexNDC(0.35, 0.75, Form("= %.2f #pm %.2f", result[0], result[1]));
        }
        //
        gROOT->ProcessLine(Form(".! mkdir -p %s", image_folder.c_str()));
        c1->SaveAs((image_folder + "/measure_breakdown_0.pdf").c_str());
        delete c1;
    }

    delete log_graph;
    return result;
}

std::array<float, 2> breakdown_voltage::measure_breakdown_1(TGraphErrors *gTarget, std::string image_folder, float threshold_guess)
{
    //  Method 1:
    //  TODO: write what it does

    //  Create result container
    std::array<float, 2> result = {-1, -1};

    //  Make the log-derivative graph
    auto log_graph = graphutils::log(gTarget);
    auto ldv_graph = graphutils::derivate(log_graph);

    //  Calculate a guess for the Vbd and have the measured baseline value
    auto guess_result = breakdown_voltage::find_vbd_guess(log_graph, threshold_guess);
    auto vbd_guess = guess_result[0];
    auto mean_baseline = guess_result[1];

    TF1 *ffull = new TF1("ffull", "[0]*TMath::Landau(x,[1],[2],0)*TMath::Gaus(x,[1],[3])");
    ffull->SetParameter(0, 0.5);
    ffull->SetParameter(1, vbd_guess);
    ffull->SetParameter(2, 1);
    ffull->SetParameter(3, 0.5);
    ldv_graph->Fit(ffull, "Q", "SAME", vbd_guess - 2, vbd_guess + 2);
    result[0] = ffull->GetParameter(1);
    result[1] = ffull->GetParError(1);

    if (!image_folder.empty())
    {
        //  Create Canvas with result
        TCanvas *c1 = new TCanvas("c1", "c1", 600, 500);
        gPad->SetMargin(0.15, 0.05, 0.10, 0.05);

        //  Set-up graph graphics
        ldv_graph->GetXaxis()->SetTitle("bias voltage (V)");
        ldv_graph->GetYaxis()->SetTitle("log of current");
        ldv_graph->SetMarkerStyle(20);
        ldv_graph->SetMarkerColor(kAzure - 3);
        ldv_graph->GetXaxis()->SetRangeUser(result[0] - 6, result[0] + 3);
        ldv_graph->GetYaxis()->SetRangeUser(ffull->Eval(result[0] - 6) - 1, ffull->Eval(result[0]) + 1);
        ldv_graph->Draw("ALPE");

        //  Plot all fits
        ffull->SetRange(result[0] - 5, result[0] + 5);
        ffull->SetLineColor(kRed);
        ffull->DrawCopy("SAME");

        //  Write on canvas results
        TLatex *l1 = new TLatex();
        l1->DrawLatexNDC(0.45, 0.90, Form("V_{bd} (V)"));
        l1->DrawLatexNDC(0.18, 0.85, Form("guess"));
        l1->DrawLatexNDC(0.35, 0.85, Form("= %.2f ", vbd_guess));
        l1->DrawLatexNDC(0.18, 0.80, Form("fit"));
        l1->DrawLatexNDC(0.35, 0.80, Form("= %.2f #pm %.2f", result[0], result[1]));
        //
        gROOT->ProcessLine(Form(".! mkdir -p %s", image_folder.c_str()));
        c1->SaveAs((image_folder + "/measure_breakdown_1.pdf").c_str());
        delete c1;
    }
    delete log_graph;
    delete ldv_graph;
    return result;
}

std::array<float, 2> breakdown_voltage::measure_breakdown_2(TGraphErrors *gTarget, std::string image_folder, float threshold_guess)
{
    //  Method 2:
    //  TODO: write what it does

    //  Create result container
    std::array<float, 2> result = {-1, -1};

    //  Make the log-derivative graph
    auto log_graph = graphutils::log(gTarget);
    auto ldv_graph = graphutils::derivate(log_graph);
    auto ild_graph = graphutils::power(ldv_graph, -1);

    //  Calculate a guess for the Vbd and have the measured baseline value
    auto guess_result = breakdown_voltage::find_vbd_guess(log_graph, threshold_guess);
    auto vbd_guess = guess_result[0];
    auto mean_baseline = guess_result[1];

    //  Fit the graph
    TF1 *pol1_after = new TF1("pol1_after", "[0]+x*[1]");
    ild_graph->Fit(pol1_after, "Q", "SAME", vbd_guess + 0.5, vbd_guess + 3.5);

    //  Get intercept
    auto intercept = get_intercept_x(pol1_after);
    result[0] = intercept[0][0];
    result[1] = intercept[0][1];

    if (!image_folder.empty())
    {
        //  Create Canvas with result
        TCanvas *c1 = new TCanvas("c1", "c1", 600, 500);
        gPad->SetMargin(0.15, 0.05, 0.10, 0.05);

        //  Set-up graph graphics
        ild_graph->GetXaxis()->SetTitle("bias voltage (V)");
        ild_graph->GetYaxis()->SetTitle("inverse log of current");
        ild_graph->SetMarkerStyle(20);
        ild_graph->SetMarkerColor(kAzure - 3);
        ild_graph->GetXaxis()->SetRangeUser(result[0] - 3, result[0] + 6);
        ild_graph->GetYaxis()->SetRangeUser(pol1_after->Eval(result[0] - 1), pol1_after->Eval(result[0] + 4) + 1);
        ild_graph->Draw("ALPE");

        //  Plot all fits
        pol1_after->SetRange(result[0] - 5, result[0] + 5);
        pol1_after->SetLineColor(kRed);
        pol1_after->DrawCopy("SAME");

        //  Write on canvas results
        TLatex *l1 = new TLatex();
        l1->DrawLatexNDC(0.45, 0.90, Form("V_{bd} (V)"));
        l1->DrawLatexNDC(0.18, 0.85, Form("guess"));
        l1->DrawLatexNDC(0.35, 0.85, Form("= %.2f ", vbd_guess));
        l1->DrawLatexNDC(0.18, 0.80, Form("fit"));
        l1->DrawLatexNDC(0.35, 0.80, Form("= %.2f #pm %.2f", result[0], result[1]));
        //
        gROOT->ProcessLine(Form(".! mkdir -p %s", image_folder.c_str()));
        c1->SaveAs((image_folder + "/measure_breakdown_2.pdf").c_str());
        delete c1;
    }

    delete log_graph;
    delete ldv_graph;
    delete ild_graph;
    return result;
}

std::array<float, 2> breakdown_voltage::measure_breakdown_3(TGraphErrors *gTarget, std::string image_folder, float threshold_guess)
{
    //  Method 3:
    //  TODO: write what it does

    //  Create result container
    std::array<float, 2> result = {-1, -1};

    //  Make the log-derivative graph
    auto log_graph = graphutils::log(gTarget);
    auto ldv_graph = graphutils::derivate(log_graph);
    auto ld2_graph = graphutils::derivate(ldv_graph);

    //  Calculate a guess for the Vbd and have the measured baseline value
    auto guess_result = breakdown_voltage::find_vbd_guess(log_graph, threshold_guess);
    auto vbd_guess = guess_result[0];
    auto mean_baseline = guess_result[1];

    //  Fit the graph
    TF1 *gaus_peak = new TF1("gaus_peak", "gaus(0)");
    ld2_graph->Fit(gaus_peak, "Q", "SAME", vbd_guess - 3.5, vbd_guess + 3.5);
    result[0] = gaus_peak->GetParameter(1);
    result[1] = gaus_peak->GetParError(1);

    if (!image_folder.empty())
    {
        //  Create Canvas with result
        TCanvas *c1 = new TCanvas("c1", "c1", 600, 500);
        gPad->SetMargin(0.15, 0.05, 0.10, 0.05);

        //  Set-up graph graphics
        ld2_graph->GetXaxis()->SetTitle("bias voltage (V)");
        ld2_graph->GetYaxis()->SetTitle("double derivative of log of current");
        ld2_graph->SetMarkerStyle(20);
        ld2_graph->SetMarkerColor(kAzure - 3);
        ld2_graph->GetXaxis()->SetRangeUser(result[0] - 4, result[0] + 4);
        ld2_graph->GetYaxis()->SetRangeUser(gaus_peak->Eval(result[0] + 4) - 1, gaus_peak->Eval(result[0]) + 1);
        ld2_graph->Draw("ALPE");

        //  Plot all fits
        gaus_peak->SetRange(result[0] - 4, result[0] + 4);
        gaus_peak->SetLineColor(kRed);
        gaus_peak->DrawCopy("SAME");

        //  Write on canvas results
        TLatex *l1 = new TLatex();
        l1->DrawLatexNDC(0.45, 0.90, Form("V_{bd} (V)"));
        l1->DrawLatexNDC(0.18, 0.85, Form("guess"));
        l1->DrawLatexNDC(0.35, 0.85, Form("= %.2f ", vbd_guess));
        l1->DrawLatexNDC(0.18, 0.80, Form("fit"));
        l1->DrawLatexNDC(0.35, 0.80, Form("= %.2f #pm %.2f", result[0], result[1]));
        //
        gROOT->ProcessLine(Form(".! mkdir -p %s", image_folder.c_str()));
        c1->SaveAs((image_folder + "/measure_breakdown_3.pdf").c_str());
        delete c1;
    }

    delete log_graph;
    delete ldv_graph;
    delete ld2_graph;
    return result;
}

std::array<float, 2> breakdown_voltage::measure_breakdown_4(TGraphErrors *gTarget, std::string image_folder, float threshold_guess)
{
    //  Method 4:
    //  TODO: write what it does

    //  Create result container
    std::array<float, 2> result = {-1, -1};

    /*
        //  Make the log graph
        auto log_graph = graphutils::log(gTarget);

        //  Calculate a guess for the Vbd and have the measured baseline value
        auto guess_result = breakdown_voltage::find_vbd_guess(log_graph, threshold_guess);
        auto vbd_guess = guess_result[0];
        auto mean_baseline = guess_result[1];

        std::array<float, 2> result = {-1, -1};
        auto log_graph = graphutils::log(gTarget);
        //  Calculate a guess for the Vbd and have the measured baseline value
        auto vbd_guess_baseline = find_vbd_guess(log_graph);
        vbd_guess = get<0>(vbd_guess_baseline);
        auto mean_baseline = get<1>(vbd_guess_baseline);
        TF1 *pol1_before = new TF1("pol1_before", "[0]+x*[1]");
        TF1 *pol1_after = new TF1("pol1_after", "[0]+x*[1]+x*x*[2]");
        pol1_before->SetParameter(0, mean_baseline);
        pol1_before->SetParameter(1, 0.);
        pol1_after->SetParameter(0, 0.);
        pol1_after->SetParameter(1, 4);
        pol1_after->SetParLimits(2, -1.e3, 0.);
        pol1_after->SetParameter(2, -0.03);
        log_graph->Fit(pol1_before, "Q", "", vbd_guess - 5, vbd_guess - 0.3);
        log_graph->Fit(pol1_after, "Q", "", vbd_guess + 0.3, vbd_guess + 2.3);
        auto intercept_before = pol1_before->GetParameter(0);
        auto angularc_before = pol1_before->GetParameter(1);
        auto c0_after = pol1_after->GetParameter(0);
        auto c1_after = pol1_after->GetParameter(1);
        auto c2_after = pol1_after->GetParameter(2);
        auto delta = sqrt((c1_after - angularc_before) * (c1_after - angularc_before) - 4 * c2_after * (c0_after - intercept_before));
        auto intercept_before_econtrib = pol1_before->GetParError(0) / (delta);
        auto angularc_before_econtrib = pol1_before->GetParError(1) * ((c1_after - angularc_before) / (delta) + 1) / (2 * c2_after);
        auto c0_after_econtrib = pol1_after->GetParError(0) / (delta);
        auto c1_after_econtrib = -pol1_after->GetParError(1) * ((c1_after - angularc_before) / (delta) + 1) / (2 * c2_after);
        auto c2_after_econtrib = pol1_after->GetParError(2) * ((c0_after - intercept_before) / (c2_after * delta) - (-delta - c1_after + angularc_before) / (2 * c2_after * c2_after));
        auto intercept_val = result[0] = (-(c1_after - angularc_before) + delta) / (2 * c2_after);
        auto intercept_err = result[1] = sqrt(intercept_before_econtrib * intercept_before_econtrib + angularc_before_econtrib * angularc_before_econtrib + c0_after_econtrib * c0_after_econtrib + c1_after_econtrib * c1_after_econtrib + c2_after_econtrib * c2_after_econtrib);
        TF1 *ffull = new TF1("ffull", "(0.5-0.5*TMath::Erf(((x-[0])/[1])))*([2]+x*[3])+(0.5+0.5*TMath::Erf(((x-[0])/[1])))*([4]+x*[5]+x*x*[6])");
        ffull->SetParameter(0, result[0] + 0.2);
        ffull->SetParLimits(0, intercept_val - 0.3, intercept_val + 5.);
        ffull->SetParameter(1, 0.2);
        ffull->SetParLimits(1, 0.0001, 0.3);
        ffull->FixParameter(2, pol1_before->GetParameter(0));
        ffull->FixParameter(3, pol1_before->GetParameter(1));
        ffull->FixParameter(4, pol1_after->GetParameter(0));
        ffull->FixParameter(5, pol1_after->GetParameter(1));
        ffull->FixParameter(6, pol1_after->GetParameter(2));
        log_graph->Fit(ffull, "Q", "SAME", result[0] - 4, result[0] + 4.);
        // ffull->SetParLimits(0, intercept_val + ffull->GetParameter(1), intercept_val + 5.);
        // log_graph->Fit(ffull, "RQ", "SAME", result[0] - 4, result[0] + 3);
        result[0] = ffull->GetParameter(0) - ffull->GetParameter(1);
        result[1] = sqrt(ffull->GetParError(0) * ffull->GetParError(0) + ffull->GetParameter(1) * ffull->GetParameter(1));
        if (graphics)
        {
            TCanvas *c1 = new TCanvas(gTarget->GetName(), "c1", 600, 500);
            log_graph->Draw("ALP");
            log_graph->GetXaxis()->SetRangeUser(result[0] - 6, result[0] + 3);
            auto max = ffull->Eval(result[0] + 2);
            max = max > 0 ? max * 1.1 : max * 0.9;
            log_graph->SetMaximum(max);
            pol1_before->SetRange(result[0] - 5, result[0] + 5);
            pol1_before->SetLineColor(kBlue);
            pol1_before->DrawCopy("SAME");
            pol1_after->SetRange(result[0] - 5, result[0] + 5);
            pol1_after->SetLineColor(kGreen);
            pol1_after->DrawCopy("SAME");
            ffull->SetRange(result[0] - 5, result[0] + 5);
            ffull->DrawCopy("SAME");
            TLatex *l1 = new TLatex();
            l1->DrawLatexNDC(0.2, 0.85, Form("V_{bd} = %.2f #pm %.2f", result[0], result[1]));
            l1->DrawLatexNDC(0.2, 0.80, Form("V_{bd} (int.) = %.2f #pm %.2f", intercept_val, intercept_err));
            //
            TString basedir = get_environment_variable("SIPM4EIC_WORKDIR");
            basedir += TString("/plots/vbdcheck/");
            gROOT->ProcessLine(Form(".! mkdir -p %s", basedir.Data()));
            c1->SaveAs(basedir + gTarget->GetName() + TString(".pdf"));
            delete c1;
        }
        delete log_graph;
        */
    return result;
}

std::array<float, 2> measure_breakdown_5(TGraphErrors *gTarget, std::string image_folder, float threshold_guess)
{
    //  Method 5:
    //  TODO: write what it does

    //  Create result container
    std::array<float, 2> result = {-1, -1};

    /*

    //  Make the log-derivative graph
    auto log_graph = graphutils::log(gTarget);
    auto ldv_graph = graphutils::derivate(log_graph);
    auto ld2_graph = graphutils::derivate(ldv_graph);

    //  Calculate a guess for the Vbd and have the measured baseline value
    auto guess_result = breakdown_voltage::find_vbd_guess(log_graph, threshold_guess);
    auto vbd_guess = guess_result[0];
    auto mean_baseline = guess_result[1];

    //  Fit the graph
    TF1 *gaus_peak = new TF1("gaus_peak", "gaus(0)");
    ld2_graph->Fit(gaus_peak, "Q", "SAME", vbd_guess - 3.5, vbd_guess + 3.5);
    result[0] = gaus_peak->GetParameter(1);
    result[1] = gaus_peak->GetParError(1);

    if (!image_folder.empty())
    {
        //  Create Canvas with result
        TCanvas *c1 = new TCanvas("c1", "c1", 600, 500);
        gPad->SetMargin(0.15, 0.05, 0.10, 0.05);

        //  Set-up graph graphics
        ld2_graph->GetXaxis()->SetTitle("bias voltage (V)");
        ld2_graph->GetYaxis()->SetTitle("double derivative of log of current");
        ld2_graph->SetMarkerStyle(20);
        ld2_graph->SetMarkerColor(kAzure - 3);
        ld2_graph->GetXaxis()->SetRangeUser(result[0] - 4, result[0] + 4);
        ld2_graph->GetYaxis()->SetRangeUser(gaus_peak->Eval(result[0] + 4) - 1, gaus_peak->Eval(result[0]) + 1);
        ld2_graph->Draw("ALPE");

        //  Plot all fits
        gaus_peak->SetRange(result[0] - 4, result[0] + 4);
        gaus_peak->SetLineColor(kRed);
        gaus_peak->DrawCopy("SAME");

        //  Write on canvas results
        TLatex *l1 = new TLatex();
        l1->DrawLatexNDC(0.45, 0.90, Form("V_{bd} (V)"));
        l1->DrawLatexNDC(0.18, 0.85, Form("guess"));
        l1->DrawLatexNDC(0.35, 0.85, Form("= %.2f ", vbd_guess));
        l1->DrawLatexNDC(0.18, 0.80, Form("fit"));
        l1->DrawLatexNDC(0.35, 0.80, Form("= %.2f #pm %.2f", result[0], result[1]));
        //
        gROOT->ProcessLine(Form(".! mkdir -p %s", image_folder.c_str()));
        c1->SaveAs((image_folder + "/measure_breakdown_3.pdf").c_str());
        delete c1;
    }

    delete log_graph;
    delete ldv_graph;
    delete ld2_graph;
    return result;

    std::array<float, 2> result = {-1, -1};
    auto log_graph = graphutils::log(gTarget);
    TF1 *pol1_before = new TF1("pol1_before", "[0]+x*[1]");
    TF1 *pol1_after = new TF1("pol1_after", "[0]*(x-[1])^([2])");
    pol1_after->SetParameter(1, vbd_guess);
    log_graph->Fit(pol1_before, "MERQ", "", vbd_guess - 5, vbd_guess - 0.2);
    log_graph->Fit(pol1_after, "MERQ", "", vbd_guess + 0.2, vbd_guess + 5);
    TF1 *ffull = new TF1("ffull", "([1]+x*[2])+(x>[0])*(+[3]*(x-[4])^([5]))");
    ffull->SetParameter(0, vbd_guess);
    ffull->FixParameter(1, pol1_before->GetParameter(0));
    ffull->FixParameter(2, pol1_before->GetParameter(1));
    ffull->FixParameter(3, pol1_after->GetParameter(0));
    ffull->FixParameter(4, pol1_after->GetParameter(1));
    ffull->FixParameter(5, pol1_after->GetParameter(2));
    log_graph->Fit(ffull, "MERQ", "SAME", vbd_guess - 3, vbd_guess + 3);
    result[0] = 48; // ffull->GetParameter(0) - ffull->GetParameter(1);
    result[1] = sqrt(ffull->GetParError(0) * ffull->GetParError(0) + ffull->GetParameter(1) * ffull->GetParameter(1));
    if (graphics)
    {
        TCanvas *c1 = new TCanvas();
        log_graph->Draw("ALP");
        log_graph->GetXaxis()->SetRangeUser(result[0] - 2, result[0] + 2);
        pol1_before->SetRange(result[0] - 5, result[0] + 5);
        pol1_before->SetLineColor(kBlue);
        pol1_before->DrawCopy("SAME");
        pol1_after->SetRange(result[0] - 5, result[0] + 5);
        pol1_after->SetLineColor(kGreen);
        pol1_after->DrawCopy("SAME");
        ffull->SetRange(result[0] - 5, result[0] + 5);
        ffull->DrawCopy("SAME");
        TLatex *l1 = new TLatex();
        l1->DrawLatexNDC(0.2, 0.85, Form("V_{bd} = %.2f #pm %.2f", result[0], result[1]));
        //
        TString basedir = get_environment_variable("SIPM4EIC_WORKDIR");
        basedir += TString("/plots/vbdcheck/");
        gROOT->ProcessLine(Form(".! mkdir -p %s", basedir.Data()));
        c1->SaveAs(basedir + gTarget->GetName() + TString(".pdf"));
        delete c1;
    }
    delete log_graph;
    result[0] = ffull->GetParameter(0);
    result[1] = sqrt(ffull->GetParError(0) * ffull->GetParError(0));
        */
    return result;
}

/*

//  List of methods to determine the breakdwon voltage in a I-V curve
//  methods taken from https://arxiv.org/abs/1606.07805
//
//  !TODO: Clean-up graph utils and make new general util repository
//
typedef std::map<std::string, std::map<std::string, std::map<std::string, std::array<double, 2>>>> breakdown_voltages_type;
typedef std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, std::array<double, 2>>>>> breakdown_voltages_sensors_type;
//
//
std::array<float, 2>
measure_breakdown_4(TGraphErrors *gTarget, float vbd_guess, bool graphics = false)
{
    std::array<float, 2> result = {-1, -1};
    auto log_graph = graphutils::log(gTarget);
    //  Calculate a guess for the Vbd and have the measured baseline value
    auto vbd_guess_baseline = find_vbd_guess(log_graph);
    vbd_guess = get<0>(vbd_guess_baseline);
    auto mean_baseline = get<1>(vbd_guess_baseline);
    TF1 *pol1_before = new TF1("pol1_before", "[0]+x*[1]");
    TF1 *pol1_after = new TF1("pol1_after", "[0]+x*[1]+x*x*[2]");
    pol1_before->SetParameter(0, mean_baseline);
    pol1_before->SetParameter(1, 0.);
    pol1_after->SetParameter(0, 0.);
    pol1_after->SetParameter(1, 4);
    pol1_after->SetParLimits(2, -1.e3, 0.);
    pol1_after->SetParameter(2, -0.03);
    log_graph->Fit(pol1_before, "Q", "", vbd_guess - 5, vbd_guess - 0.3);
    log_graph->Fit(pol1_after, "Q", "", vbd_guess + 0.3, vbd_guess + 2.3);
    auto intercept_before = pol1_before->GetParameter(0);
    auto angularc_before = pol1_before->GetParameter(1);
    auto c0_after = pol1_after->GetParameter(0);
    auto c1_after = pol1_after->GetParameter(1);
    auto c2_after = pol1_after->GetParameter(2);
    auto delta = sqrt((c1_after - angularc_before) * (c1_after - angularc_before) - 4 * c2_after * (c0_after - intercept_before));
    auto intercept_before_econtrib = pol1_before->GetParError(0) / (delta);
    auto angularc_before_econtrib = pol1_before->GetParError(1) * ((c1_after - angularc_before) / (delta) + 1) / (2 * c2_after);
    auto c0_after_econtrib = pol1_after->GetParError(0) / (delta);
    auto c1_after_econtrib = -pol1_after->GetParError(1) * ((c1_after - angularc_before) / (delta) + 1) / (2 * c2_after);
    auto c2_after_econtrib = pol1_after->GetParError(2) * ((c0_after - intercept_before) / (c2_after * delta) - (-delta - c1_after + angularc_before) / (2 * c2_after * c2_after));
    auto intercept_val = result[0] = (-(c1_after - angularc_before) + delta) / (2 * c2_after);
    auto intercept_err = result[1] = sqrt(intercept_before_econtrib * intercept_before_econtrib + angularc_before_econtrib * angularc_before_econtrib + c0_after_econtrib * c0_after_econtrib + c1_after_econtrib * c1_after_econtrib + c2_after_econtrib * c2_after_econtrib);
    TF1 *ffull = new TF1("ffull", "(0.5-0.5*TMath::Erf(((x-[0])/[1])))*([2]+x*[3])+(0.5+0.5*TMath::Erf(((x-[0])/[1])))*([4]+x*[5]+x*x*[6])");
    ffull->SetParameter(0, result[0] + 0.2);
    ffull->SetParLimits(0, intercept_val - 0.3, intercept_val + 5.);
    ffull->SetParameter(1, 0.2);
    ffull->SetParLimits(1, 0.0001, 0.3);
    ffull->FixParameter(2, pol1_before->GetParameter(0));
    ffull->FixParameter(3, pol1_before->GetParameter(1));
    ffull->FixParameter(4, pol1_after->GetParameter(0));
    ffull->FixParameter(5, pol1_after->GetParameter(1));
    ffull->FixParameter(6, pol1_after->GetParameter(2));
    log_graph->Fit(ffull, "Q", "SAME", result[0] - 4, result[0] + 4.);
    // ffull->SetParLimits(0, intercept_val + ffull->GetParameter(1), intercept_val + 5.);
    // log_graph->Fit(ffull, "RQ", "SAME", result[0] - 4, result[0] + 3);
    result[0] = ffull->GetParameter(0) - ffull->GetParameter(1);
    result[1] = sqrt(ffull->GetParError(0) * ffull->GetParError(0) + ffull->GetParameter(1) * ffull->GetParameter(1));
    if (graphics)
    {
        TCanvas *c1 = new TCanvas(gTarget->GetName(), "c1", 600, 500);
        log_graph->Draw("ALP");
        log_graph->GetXaxis()->SetRangeUser(result[0] - 6, result[0] + 3);
        auto max = ffull->Eval(result[0] + 2);
        max = max > 0 ? max * 1.1 : max * 0.9;
        log_graph->SetMaximum(max);
        pol1_before->SetRange(result[0] - 5, result[0] + 5);
        pol1_before->SetLineColor(kBlue);
        pol1_before->DrawCopy("SAME");
        pol1_after->SetRange(result[0] - 5, result[0] + 5);
        pol1_after->SetLineColor(kGreen);
        pol1_after->DrawCopy("SAME");
        ffull->SetRange(result[0] - 5, result[0] + 5);
        ffull->DrawCopy("SAME");
        TLatex *l1 = new TLatex();
        l1->DrawLatexNDC(0.2, 0.85, Form("V_{bd} = %.2f #pm %.2f", result[0], result[1]));
        l1->DrawLatexNDC(0.2, 0.80, Form("V_{bd} (int.) = %.2f #pm %.2f", intercept_val, intercept_err));
        //
        TString basedir = get_environment_variable("SIPM4EIC_WORKDIR");
        basedir += TString("/plots/vbdcheck/");
        gROOT->ProcessLine(Form(".! mkdir -p %s", basedir.Data()));
        c1->SaveAs(basedir + gTarget->GetName() + TString(".pdf"));
        delete c1;
    }
    delete log_graph;
    return result;
}
//
//
std::array<float, 2>
measure_breakdown_6(TGraphErrors *gTarget, float vbd_guess, bool graphics = false)
{
    std::array<float, 2> result = {-1, -1};
    auto current_graph = (TGraphErrors *)gTarget->Clone("tmp");
    TF1 *fprebefore = new TF1("fprebefore", "[0]");
    TF1 *pol1_before = new TF1("pol1_before", "[0]+x*[1]");
    TF1 *pol1_after = new TF1("pol1_after", "(TMath::Sign([3],(x-[0])))*(TMath::Power(TMath::Abs(x-[0])/([1]),[2]))");
    pol1_before->FixParameter(0, fprebefore->GetParameter(0));
    pol1_after->FixParameter(0, vbd_guess);
    pol1_after->SetParLimits(1, 1.e-12, 1.e1);
    pol1_after->SetParameter(1, .7);
    pol1_after->SetParLimits(2, 1.e-12, 1.e1);
    pol1_after->SetParameter(2, .1);
    pol1_after->SetParameter(3, fprebefore->GetParameter(0));
    pol1_after->FixParameter(0, vbd_guess);
    pol1_after->SetParameter(1, 7.17717e-01);
    pol1_after->SetParameter(2, 1.);
    pol1_after->SetParameter(3, 4.59948e-09);
    current_graph->Fit(pol1_before, "MERQ", "", vbd_guess - 5, vbd_guess - 0.2);
    current_graph->Fit(pol1_after, "MERQ", "", vbd_guess + 0.2, vbd_guess + 5);
    TF1 *ffull = new TF1("ffull", "(0.5-0.5*TMath::Erf((x-[0])/[1]))*([2]+x*[3])+(0.5+0.5*TMath::Erf(((x-[0])/[1])))*([2]+x*[3]+(TMath::Sign([7],(x-[4])))*(TMath::Power(TMath::Abs(x-[4])/([5]),[6])))");
    ffull->SetParameter(0, vbd_guess);
    ffull->SetParameter(1, 0.2);
    ffull->SetParLimits(1, 0.0001, 0.3);
    ffull->FixParameter(2, pol1_before->GetParameter(0));
    ffull->FixParameter(3, pol1_before->GetParameter(1));
    ffull->FixParameter(4, pol1_after->GetParameter(0));
    ffull->FixParameter(5, pol1_after->GetParameter(1));
    ffull->FixParameter(6, pol1_after->GetParameter(2));
    ffull->FixParameter(7, pol1_after->GetParameter(3));
    current_graph->Fit(ffull, "MERQ", "SAME", vbd_guess - 5, vbd_guess + 5);
    result[0] = ffull->GetParameter(0) - ffull->GetParameter(1);
    result[1] = sqrt(ffull->GetParError(0) * ffull->GetParError(0) + ffull->GetParameter(1) * ffull->GetParameter(1));
    if (graphics)
    {
        TCanvas *c1 = new TCanvas();
        gPad->SetLogy();
        current_graph->Draw("ALP");
        current_graph->GetXaxis()->SetRangeUser(result[0] - 5, result[0] + 5);
        pol1_before->SetRange(result[0] - 5, result[0] + 5);
        pol1_before->SetLineColor(kBlue);
        pol1_before->DrawCopy("SAME");
        pol1_after->SetRange(result[0] - 5, result[0] + 5);
        pol1_after->SetLineColor(kGreen);
        pol1_after->DrawCopy("SAME");
        ffull->SetRange(result[0] - 5, result[0] + 5);
        ffull->DrawCopy("SAME");
        TLatex *l1 = new TLatex();
        l1->DrawLatexNDC(0.2, 0.85, Form("V_{bd} = %.2f #pm %.2f", result[0], result[1]));
        //
        TString basedir = get_environment_variable("SIPM4EIC_WORKDIR");
        basedir += TString("/plots/vbdcheck/");
        gROOT->ProcessLine(Form(".! mkdir -p %s", basedir.Data()));
        c1->SaveAs(basedir + gTarget->GetName() + TString(".pdf"));
        delete c1;
    }
    delete current_graph;
    result[0] = ffull->GetParameter(0);
    result[1] = sqrt(ffull->GetParError(0) * ffull->GetParError(0));
    return result;
}
//
std::array<float, 2>
measure_breakdown_100(TGraphErrors *gTarget, float vbd_guess, bool graphics = false)
{
    // Generate result pair
    std::array<float, 2> result = {-1, -1};
    // Geenrate local copy of graph
    auto current_graph = (TGraphErrors *)gTarget->Clone();
    // Linear fit to find with extrapolation the crossing of X-axis
    TF1 *pol1_after = new TF1("pol1_after", "(x-[0])/[1]");
    TF1 *ffull = new TF1("ffull", "(x-[0])/[1]+[2]/(x-[3])");
    ffull->SetParLimits(2, 0., 1.e8);
    // Function prepping
    // After fitting
    pol1_after->FixParameter(0, vbd_guess);
    pol1_after->SetParameter(1, 1. / (current_graph->Eval(vbd_guess + 3) - current_graph->Eval(vbd_guess + 2)));
    current_graph->Fit(pol1_after, "MEQ", "", vbd_guess + 2.0, vbd_guess + 2.5);
    pol1_after->SetParLimits(0, vbd_guess - 3, vbd_guess + 3);
    current_graph->Fit(pol1_after, "MEQ", "", vbd_guess + 2.0, vbd_guess + 2.5);
    // Full fitting
    ffull->FixParameter(0, pol1_after->GetParameter(0));
    ffull->FixParameter(1, pol1_after->GetParameter(1));
    ffull->SetParameter(2, 1.e+05);
    ffull->FixParameter(3, pol1_after->GetParameter(0) + 0.5);
    current_graph->Fit(ffull, "MEQ", "", vbd_guess + 0.5, vbd_guess + 2.5);
    ffull->SetParLimits(0, vbd_guess - 3, vbd_guess + 3);
    ffull->ReleaseParameter(1);
    ffull->SetParLimits(3, vbd_guess - 3, vbd_guess + 3);
    current_graph->Fit(ffull, "MEQ", "", vbd_guess + 0.5, vbd_guess + 2.5);
    // Set Parameters for show
    pol1_after->SetParameter(0, ffull->GetParameter(0));
    pol1_after->SetParameter(1, ffull->GetParameter(1));
    // Save result
    result[0] = ffull->GetParameter(0);
    result[1] = ffull->GetParError(0);
    // Plot if requested
    if (graphics)
    {
        TCanvas *c1 = new TCanvas();
        gPad->SetLogy();
        current_graph->Draw("ALP");
        current_graph->GetXaxis()->SetRangeUser(result[0] - 5, result[0] + 5);
        pol1_after->SetRange(result[0] - 5, result[0] + 5);
        pol1_after->SetLineColor(kBlue);
        pol1_after->DrawCopy("SAME");
        ffull->SetRange(result[0] - 5, result[0] + 5);
        ffull->SetLineColor(kRed);
        ffull->DrawCopy("SAME");
        // Declare found Vbd value
        TLatex *l1 = new TLatex();
        l1->DrawLatexNDC(0.2, 0.85, Form("V_{bd} = %.2f #pm %.2f", result[0], result[1]));
        // Save plot
        TString basedir = get_environment_variable("SIPM4EIC_WORKDIR");
        basedir += TString("/plots/vbdcheck/");
        gROOT->ProcessLine(Form(".! mkdir -p %s", basedir.Data()));
        c1->SaveAs(basedir + gTarget->GetName() + TString(".pdf"));
        delete c1;
    }
    delete current_graph;
    return result;
}
//
template <int mode_select = 6>
std::array<float, 2>
measure_breakdown(TGraphErrors *gTarget, float vbd_guess, bool graphics = true)
{
    switch (mode_select)
    {
        // Breakdown Voltage from IV curves
    case 0:
        return measure_breakdown_0(gTarget, vbd_guess, graphics);
        break;
    case 1:
        return measure_breakdown_1(gTarget, vbd_guess, graphics);
        break;
    case 2:
        return measure_breakdown_2(gTarget, vbd_guess, graphics);
        break;
    case 3:
        return measure_breakdown_3(gTarget, vbd_guess, graphics);
        break;
    case 4:
        return measure_breakdown_4(gTarget, vbd_guess, graphics);
        break;
    case 5:
        return measure_breakdown_5(gTarget, vbd_guess, graphics);
        break;
    case 6:
        return measure_breakdown_6(gTarget, vbd_guess, graphics);
        break;
        // Breakdown Voltage from DCR curves
    case 100:
        return measure_breakdown_100(gTarget, vbd_guess, graphics);
        break;
    default:
        return measure_breakdown_4(gTarget, vbd_guess, graphics);
        break;
    }
}
//
template <int mode_select = 1000>
std::array<float, 2>
measure_breakdown(TGraphErrors *gTarget, float vbd_guess, int ncycles, float tolerance = 1.e-2)
{
    auto initial_guess = measure_breakdown<mode_select>(gTarget, vbd_guess);
    for (int iTer = 0; iTer < ncycles; iTer++)
    {
        auto new_guess = measure_breakdown<mode_select>(gTarget, get<0>(initial_guess));
        if (fabs(get<0>(new_guess) - get<0>(initial_guess)) < tolerance)
        {
            return new_guess;
        }
        initial_guess = new_guess;
    }
    return initial_guess;
}
//
template <bool recalulate_vbd = true,
          int method = 4>
TGraphErrors *
shift_with_vbd(TGraphErrors *gTarget, float vbd_guess)
{
    auto new_target = (TGraphErrors *)gTarget->Clone(gTarget->GetName() + TString("_tmp"));
    auto vbd_shift_by = 0.;
    if (recalulate_vbd)
    {
        vbd_shift_by = get<0>(measure_breakdown<method>(new_target, vbd_guess));
    }
    else
    {
        vbd_shift_by = vbd_guess;
    }
    graphutils::x_shift(new_target, vbd_shift_by);
    return new_target;
}
//
template <bool recalulate_vbd = true>
std::array<float, 2>
overvoltage_ratio(TGraphErrors *gNumerator, TGraphErrors *gDenominator, float guess_vbd_num, float guess_vbd_den, float _xcoordinate)
{
    std::array<float, 2> result;
    std::array<float, 2> numerator;
    std::array<float, 2> denominator;
    if (recalulate_vbd)
    {
        gNumerator->SetName("gNumerator");
        gDenominator->SetName("gDenominator");
        numerator = eval_with_errors(shift_with_vbd(gNumerator, guess_vbd_num), _xcoordinate);
        denominator = eval_with_errors(shift_with_vbd(gDenominator, guess_vbd_den), _xcoordinate);
    }
    else
    {
        auto new_numerator = (TGraphErrors *)gNumerator->Clone("tmp");
        auto new_denominator = (TGraphErrors *)gDenominator->Clone("tmp");
        graphutils::x_shift(new_numerator, guess_vbd_num);
        graphutils::x_shift(new_denominator, guess_vbd_den);
        numerator = eval_with_errors(new_numerator, _xcoordinate);
        denominator = eval_with_errors(new_denominator, _xcoordinate);
    }
    result[0] = get<0>(numerator) / get<0>(denominator);
    result[1] = result[0] * ((get<1>(numerator) / get<0>(numerator)) * (get<1>(numerator) / get<0>(numerator)) + (get<1>(denominator) / get<0>(denominator)) * (get<1>(denominator) / get<0>(denominator)));
    return result;
}
//
template <bool recalulate_vbd = true>
TGraphErrors *
overvoltage_ratio(TGraphErrors *gNumerator, TGraphErrors *gDenominator, float guess_vbd_num, float guess_vbd_den, std::vector<float> _xcoordinates)
{
    TGraphErrors *result = new TGraphErrors();
    std::array<float, 2> current_ratio;
    std::array<float, 2> numerator;
    std::array<float, 2> denominator;
    auto new_numerator = (TGraphErrors *)gNumerator->Clone("gNumerator");
    auto new_denominator = (TGraphErrors *)gDenominator->Clone("gDenominator");
    if (recalulate_vbd)
    {
        shift_with_vbd(gNumerator, guess_vbd_num);
        shift_with_vbd(gDenominator, guess_vbd_den);
    }
    else
    {
        graphutils::x_shift(new_numerator, guess_vbd_num);
        graphutils::x_shift(new_denominator, guess_vbd_den);
    }
    for (auto _xcoordinate : _xcoordinates)
    {
        numerator = eval_with_errors(new_numerator, _xcoordinate);
        denominator = eval_with_errors(new_denominator, _xcoordinate);
        if (get<0>(denominator) == 0)
        {
            continue;
        }
        if ((get<0>(numerator) == get<1>(numerator)) && (get<0>(numerator) == -1))
        {
            continue;
        }
        if ((get<0>(denominator) == get<1>(denominator)) && (get<0>(denominator) == -1))
        {
            continue;
        }
        get<0>(current_ratio) = get<0>(numerator) / get<0>(denominator);
        auto numerator_err = (get<1>(numerator) / get<0>(numerator));
        auto denominator_err = (get<1>(denominator) / get<0>(denominator));
        get<1>(current_ratio) = get<0>(current_ratio) * sqrt(numerator_err * numerator_err + denominator_err * denominator_err);
        auto current_point = result->GetN();
        result->SetPoint(current_point, _xcoordinate, get<0>(current_ratio));
        result->SetPointError(current_point, 0, get<1>(current_ratio));
    }
    return result;
}
//
template <bool recalulate_vbd = true>
TGraphErrors *
overvoltage_difference(TGraphErrors *gNumerator, TGraphErrors *gDenominator, float guess_vbd_num, float guess_vbd_den, std::vector<float> _xcoordinates)
{
    TGraphErrors *result = new TGraphErrors();
    std::array<float, 2> current_ratio;
    std::array<float, 2> numerator;
    std::array<float, 2> denominator;
    auto new_numerator = (TGraphErrors *)gNumerator->Clone("gNumerator");
    auto new_denominator = (TGraphErrors *)gDenominator->Clone("gDenominator");
    if (recalulate_vbd)
    {
        shift_with_vbd(gNumerator, guess_vbd_num);
        shift_with_vbd(gDenominator, guess_vbd_den);
    }
    else
    {
        graphutils::x_shift(new_numerator, guess_vbd_num);
        graphutils::x_shift(new_denominator, guess_vbd_den);
    }
    for (auto _xcoordinate : _xcoordinates)
    {
        numerator = eval_with_errors(new_numerator, _xcoordinate);
        denominator = eval_with_errors(new_denominator, _xcoordinate);
        if ((get<0>(numerator) == get<1>(numerator)) && (get<0>(numerator) == -1))
        {
            continue;
        }
        if ((get<0>(denominator) == get<1>(denominator)) && (get<0>(denominator) == -1))
        {
            continue;
        }
        get<0>(current_ratio) = get<0>(numerator) - get<0>(denominator);
        get<1>(current_ratio) = sqrt(get<1>(numerator) * get<1>(numerator) + get<1>(numerator) * get<1>(numerator));
        auto current_point = result->GetN();
        result->SetPoint(current_point, _xcoordinate, get<0>(current_ratio));
        result->SetPointError(current_point, 0, get<1>(current_ratio));
    }
    return result;
}
//
template <int method = 6>
TGraphErrors *
measure_gain(TGraphErrors *iv_graph, TGraphErrors *dcr_graph, double vbd_guess)
{
    // Make tmp
    auto current_iv = (TGraphErrors *)iv_graph->Clone("tmp_IV");
    auto current_dcr = (TGraphErrors *)dcr_graph->Clone("tmp_DCR");

    // Shift with Vbd
    auto current_vbd = method < 0 ? vbd_guess : get<0>(measure_breakdown<method>(current_iv, vbd_guess));
    graphutils::x_shift(current_iv, current_vbd);
    graphutils::x_shift(current_dcr, current_vbd);

    // Subtract surface current from IV
    float surface_current_average = 0.;
    for (int iPoint = 0; iPoint < 5; ++iPoint)
    {
        surface_current_average += current_iv->GetY()[iPoint];
    }
    graphutils::y_shift(current_iv, surface_current_average / 5.);

    // Calculate gain
    auto current_gain = graphutils::ratio(current_iv, current_dcr);
    graphutils::y_scale(current_gain, 1. / TMath::Qe());

    // Delete tmp
    delete current_iv;
    delete current_dcr;

    // return
    return current_gain;
}
//
template <int method = 6>
std::vector<TGraphErrors *>
measure_gain(std::vector<TGraphErrors *> iv_graphs, std::vector<TGraphErrors *> dcr_graphs, std::vector<double> vbd_guess)
{
    std::vector<TGraphErrors *> result;
    auto iTer = -1;
    for (auto current_iv : iv_graphs)
    {
        iTer++;
        result.push_back(measure_gain<method>(iv_graphs.at(iTer), dcr_graphs.at(iTer), vbd_guess.at(iTer)));
    }
    //
    return result;
}
//
template <int method = 6>
std::vector<TGraphErrors *>
measure_gain(std::vector<TGraphErrors *> iv_graphs, std::vector<TGraphErrors *> dcr_graphs, double vbd_guess)
{
    std::vector<TGraphErrors *> result;
    auto iTer = -1;
    for (auto current_iv : iv_graphs)
    {
        iTer++;
        result.push_back(measure_gain<method>(iv_graphs.at(iTer), dcr_graphs.at(iTer), vbd_guess));
    }
    //
    return result;
}
//
template <int method = 6>
std::vector<TGraphErrors *>
measure_gain(std::vector<TGraphErrors *> iv_graphs, std::vector<TGraphErrors *> dcr_graphs, std::map<TString, std::array<double, 2>> breakdown_voltages)
{
    std::vector<TGraphErrors *> result;
    auto iTer = -1;
    for (auto current_iv : iv_graphs)
    {
        iTer++;
        result.push_back(measure_gain<-1>(iv_graphs.at(iTer), dcr_graphs.at(iTer), get<0>(breakdown_voltages[iv_graphs.at(iTer)->GetName()])));
    }
    //
    return result;
}
//
namespace database
{
    // https://docs.google.com/document/d/1Z3A6zrad5A-7HcmOVEIBOqT7nCkWAe4-ilvuk-I16rs/edit#heading=h.a7bxpd1bw3ht
    //
    // --- Breakdown voltage
    // Storage for on-the-fly Vbd
    //  --- board - channel - step
    breakdown_voltages_type breakdown_voltages;
    //  --- sensor - step - board - channel
    breakdown_voltages_sensors_type breakdown_voltages_sensors;
    //
    std::map<std::string, double> sensor_vbd_guess = {
        {"Hamamatsu S13361-3050", 48.3},
        {"Hamamatsu S13361-3075", 47.8},
        {"Hamamatsu S14161-3050", 36.5}};
    //
    template <int method = -1>
    void calculate_breakdown_voltages()
    {
        for (auto current_status : database::all_statuses)
        {
            for (auto current_board : database::all_boards)
            {
                for (auto current_channel : database::all_channels)
                {
                    auto vbd_guess = 0.;
                    if (current_channel.find("A") != std::string::npos)
                    {
                        vbd_guess = 48.3;
                    }
                    if (current_channel.find("B") != std::string::npos)
                    {
                        vbd_guess = 47.8;
                    }
                    if (current_channel.find("C") != std::string::npos)
                    {
                        vbd_guess = 36.5;
                    }
                    // Recover IV
                    auto current_iv = database::get_iv_scan(current_board, current_channel, current_status);
                    if (!current_iv)
                        continue;
                    current_iv->SetName(TString(current_board) + "/" + TString(current_status) + "/" + TString(current_channel));
                    TString basedir = (TString(get_environment_variable("SIPM4EIC_WORKDIR")) + TString("/plots/vbdcheck/") + TString(current_board) + "/" + TString(current_status) + "/");
                    gROOT->ProcessLine(Form(".! mkdir -p %s", basedir.Data()));
                    if (method < 0)
                    {
                        if ((current_channel.find("A") != std::string::npos) && (current_status.find("NEW") != std::string::npos))
                        {
                            breakdown_voltages[current_board][current_channel][current_status] = measure_breakdown<0>(current_iv, vbd_guess);
                            breakdown_voltages_sensors[channel_to_sensor[current_channel]][current_status][current_board][current_channel] = measure_breakdown<0>(current_iv, vbd_guess);
                        }
                        else
                        {

                            breakdown_voltages[current_board][current_channel][current_status] = measure_breakdown<4>(current_iv, vbd_guess);
                            breakdown_voltages_sensors[channel_to_sensor[current_channel]][current_status][current_board][current_channel] = measure_breakdown<4>(current_iv, vbd_guess);
                        }
                    }
                    else
                    {
                        breakdown_voltages[current_board][current_channel][current_status] = measure_breakdown<method>(current_iv, vbd_guess);
                        breakdown_voltages_sensors[channel_to_sensor[current_channel]][current_status][current_board][current_channel] = measure_breakdown<method>(current_iv, vbd_guess);
                    }
                }
            }
        }
    }
    std::array<double, 2> get_breakdown_voltage(std::string board, std::string channel, std::string step)
    {
        return breakdown_voltages[board][channel][step];
    }
    breakdown_voltages_type get_all_breakdown_voltages()
    {
        return breakdown_voltages;
    }
}*/