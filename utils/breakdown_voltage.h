#pragma once
//  List of methods to determine the breakdwon voltage in a I-V curve
//  methods taken from https://arxiv.org/abs/1606.07805
//
//  !TODO: Clean-up graph utils and make new general util repository
//
typedef std::map<std::string, std::map<std::string, std::map<std::string, std::array<double, 2>>>> breakdown_voltages_type;
typedef std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, std::array<double, 2>>>>> breakdown_voltages_sensors_type;
//
std::array<double, 2> find_vbd_guess(TGraphErrors *gLogTarget)
{
    auto vbd_guess = 0.;
    auto mean_baseline = gLogTarget->GetY()[0];
    auto mean_contributors = 1;
    for (int iPnt = 1; iPnt < gLogTarget->GetN(); iPnt++)
    {
        auto current_Y = gLogTarget->GetY()[iPnt];
        if (fabs(mean_baseline - current_Y) > 0.35)
        {
            vbd_guess = gLogTarget->GetX()[iPnt-1];
            break;
        }
        mean_baseline *= mean_contributors;
        mean_baseline += current_Y;
        mean_contributors++;
        mean_baseline /= mean_contributors;
    }
    return {vbd_guess, mean_baseline};
}
//
std::pair<float, float> measure_breakdown_0(TGraphErrors *gTarget, float vbd_guess, bool graphics = false)
{
    //  Create the result pair
    std::pair<float, float> result = {-1, -1};
    //  Make the log graph
    auto log_graph = graphutils::log(gTarget);
    //  Calculate a guess for the Vbd and have the measured baseline value
    auto vbd_guess_baseline = find_vbd_guess(log_graph);
    vbd_guess = get<0>(vbd_guess_baseline);
    auto mean_baseline = get<1>(vbd_guess_baseline);
    //  Define the functions to find the Vbd
    TF1 *fbefore = new TF1("fbefore", "[0]+x*[1]");
    TF1 *fafter = new TF1("fafter", "[0]+x*[1]");
    //  Fit before and after the Vbd guess
    fbefore->SetParameter(0, mean_baseline);
    fbefore->SetParameter(1, 0.);
    log_graph->Fit(fbefore, "RQ", "", vbd_guess - 5, vbd_guess - 0.2);
    fafter->SetParameter(0, -100);
    fafter->SetParameter(1, 4);
    log_graph->Fit(fafter, "RQ", "", vbd_guess + 0.2, vbd_guess + 1.);
    auto intercept_before = fbefore->GetParameter(0);
    auto angularc_before = fbefore->GetParameter(1);
    auto intercept_after = fafter->GetParameter(0);
    auto angularc_after = fafter->GetParameter(1);
    auto numerator_contribe = sqrt(fbefore->GetParError(0) * fbefore->GetParError(0) + fafter->GetParError(0) * fafter->GetParError(0)) / (intercept_after - intercept_before);
    auto denominator_contribe = sqrt(fbefore->GetParError(1) * fbefore->GetParError(1) + fafter->GetParError(1) * fafter->GetParError(1)) / (angularc_before - angularc_after);
    auto intercept_val = get<0>(result) = (intercept_after - intercept_before) / (angularc_before - angularc_after);
    auto intercept_err = get<1>(result) = get<0>(result) * sqrt(numerator_contribe * numerator_contribe + denominator_contribe * denominator_contribe);
    TF1 *ffull = new TF1("ffull", "((0.5-0.5*TMath::Erf(((x-[0])/[1])))*([2]+x*[3])+(0.5+0.5*TMath::Erf(((x-[0])/[1])))*([4]+x*[5]))");
    ffull->SetParameter(0, intercept_val);
    ffull->SetParLimits(0, intercept_val - 0.2, intercept_val + 0.2);
    ffull->SetParameter(1, 0.0000015);
    ffull->SetParLimits(1, 0.000001, 0.01);
    ffull->FixParameter(2, fbefore->GetParameter(0));
    ffull->FixParameter(3, fbefore->GetParameter(1));
    ffull->FixParameter(4, fafter->GetParameter(0));
    ffull->FixParameter(5, fafter->GetParameter(1));
    log_graph->Fit(ffull, "RQ", "SAME", get<0>(result) - 4, get<0>(result) + 3.);
    ffull->SetParLimits(0, intercept_val + ffull->GetParameter(1), intercept_val + 0.2);
    log_graph->Fit(ffull, "RQ", "SAME", get<0>(result) - 4, get<0>(result) + +3.);
    get<0>(result) = ffull->GetParameter(0) - ffull->GetParameter(1);
    get<1>(result) = sqrt(ffull->GetParError(0) * ffull->GetParError(0) + ffull->GetParameter(1) * ffull->GetParameter(1));
    if (graphics)
    {
        TCanvas *c1 = new TCanvas(gTarget->GetName(), "c1", 600, 500);
        log_graph->Draw("ALP");
        log_graph->GetXaxis()->SetRangeUser(get<0>(result) - 6, get<0>(result) + 3);
        auto max = ffull->Eval(get<0>(result) + 2);
        max = max > 0 ? max * 1.1 : max * 0.9;
        log_graph->SetMaximum(max);
        fbefore->SetRange(get<0>(result) - 5, get<0>(result) + 5);
        fbefore->SetLineColor(kBlue);
        fbefore->DrawCopy("SAME");
        fafter->SetRange(get<0>(result) - 5, get<0>(result) + 5);
        fafter->SetLineColor(kGreen);
        fafter->DrawCopy("SAME");
        ffull->SetRange(get<0>(result) - 5, get<0>(result) + 5);
        ffull->DrawCopy("SAME");
        TLatex *l1 = new TLatex();
        l1->DrawLatexNDC(0.2, 0.85, Form("V_{bd} guess = %.2f ", vbd_guess));
        l1->DrawLatexNDC(0.2, 0.80, Form("V_{bd} = %.2f #pm %.2f", get<0>(result), get<1>(result)));
        l1->DrawLatexNDC(0.2, 0.75, Form("V_{bd} (int.) = %.2f #pm %.2f", intercept_val, intercept_err));
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
std::pair<float, float>
measure_breakdown_1(TGraphErrors *gTarget, float vbd_guess, bool graphics = true)
{
    std::pair<float, float> result = {-1, -1};
    auto log_graph = graphutils::log(gTarget);
    auto deriv_graph = graphutils::derivate(log_graph);
    auto baseline_reference = 0.;
    for (auto iPnt = 0; iPnt < 5; iPnt++)
        baseline_reference += deriv_graph->GetY()[iPnt] / 5.;
    auto vbd_reference = 0.;
    auto peak_reference = 0.;
    for (auto iPnt = 5; iPnt < 25; iPnt++)
    {
        if (peak_reference < deriv_graph->GetY()[iPnt])
        {
            vbd_reference = deriv_graph->GetX()[iPnt];
            peak_reference = deriv_graph->GetY()[iPnt];
        }
    }
    // TF1 *fpeak = new TF1("fpeak", "(([0]*[1])/(2))*TMath::Exp(([1]/2)*(2*[2]+[1]*[3]*[3]-2*x))*(1-TMath::Erf(([2]+[1]*[3]*[3]-x)/(TMath::Sqrt(2)*[3])))");
    TF1 *fpeak = new TF1("fpeak", "[0]*TMath::Landau(x,[1],[2],1)");
    TF1 *ffull = new TF1("ffull", "[0]*TMath::Landau(x,[1],[2],0)*TMath::Gaus(x,[1],[3])");
    fpeak->SetParameter(0, 1.);
    fpeak->SetParameter(1, vbd_reference);
    fpeak->SetParameter(2, .1);
    deriv_graph->Fit(fpeak, "Q", "SAME", vbd_reference - 1.5, vbd_reference + 1.5);
    deriv_graph->Fit(fpeak, "Q", "SAME", vbd_reference - 1.5, vbd_reference + 1.5);
    deriv_graph->Fit(fpeak, "Q", "SAME", vbd_reference - 1.5, vbd_reference + 1.5);
    get<0>(result) = fpeak->GetParameter(1);
    get<1>(result) = fpeak->GetParError(1);
    if (graphics)
    {
        TCanvas *c1 = new TCanvas(gTarget->GetName(), "c1", 600, 500);
        deriv_graph->Draw("ALP");
        deriv_graph->GetXaxis()->SetRangeUser(get<0>(result) - 6, get<0>(result) + 3);
        fpeak->SetRange(get<0>(result) - 6, get<0>(result) + 3);
        fpeak->SetLineColor(kBlue);
        fpeak->SetLineStyle(kDashed);
        fpeak->DrawCopy("SAME");
        ffull->SetRange(get<0>(result) - 6, get<0>(result) + 3);
        ffull->DrawCopy("SAME");
        TLatex *l1 = new TLatex();
        l1->DrawLatexNDC(0.2, 0.85, Form("V_{bd} = %.2f #pm %.2f", get<0>(result), get<1>(result)));
        l1->DrawLatexNDC(0.2, 0.80, Form("V_{bd} = %.2f", vbd_reference));
        //
        TString basedir = get_environment_variable("SIPM4EIC_WORKDIR");
        basedir += TString("/plots/vbdcheck/");
        gROOT->ProcessLine(Form(".! mkdir -p %s", basedir.Data()));
        c1->SaveAs(basedir + gTarget->GetName() + TString(".pdf"));
        delete c1;
    }
    delete log_graph;
    delete deriv_graph;
    return result;
}
//
std::pair<float, float>
measure_breakdown_2(TGraphErrors *gTarget, float vbd_guess, bool graphics = false)
{
    std::pair<float, float> result = {-1, -1};
    auto log_graph = graphutils::log(gTarget);
    auto deriv_graph = graphutils::derivate(log_graph);
    auto inverse_graph = graphutils::power(deriv_graph, -1);
    for (auto iPnt = 3; iPnt < inverse_graph->GetN() - 3; iPnt++)
    {
        auto current_x = inverse_graph->GetX()[iPnt];
        auto current_y = inverse_graph->GetY()[iPnt];
        auto prev1_y = inverse_graph->GetY()[iPnt - 1];
        auto prev2_y = inverse_graph->GetY()[iPnt - 2];
        auto prev3_y = inverse_graph->GetY()[iPnt - 2];
        auto next1_y = inverse_graph->GetY()[iPnt + 1];
        auto next2_y = inverse_graph->GetY()[iPnt + 2];
        auto next3_y = inverse_graph->GetY()[iPnt + 2];
        if (!((current_y < prev1_y) && (prev1_y < prev2_y) && (prev2_y <= prev3_y)))
            continue;
        if (!((current_y < next1_y) && (next1_y < next2_y) && (next2_y <= next3_y)))
            continue;
        vbd_guess = current_x;
        break;
    }
    TF1 *fafter = new TF1("fafter", "[0]+x*[1]");
    inverse_graph->Fit(fafter, "Q", "SAME", vbd_guess, vbd_guess + 2);
    auto qcoeff = fafter->GetParameter(0);
    auto mcoeff = fafter->GetParameter(1);
    auto qcoefe = fafter->GetParError(0);
    auto mcoefe = fafter->GetParError(1);
    get<0>(result) = -(qcoeff) / (mcoeff);
    get<1>(result) = -(qcoeff) / (mcoeff)*sqrt(((qcoefe / qcoeff) * (qcoefe / qcoeff)) + ((mcoefe / mcoeff) * (mcoefe / mcoeff)));
    if (graphics)
    {
        TCanvas *c1 = new TCanvas(gTarget->GetName(), "c1", 600, 500);
        inverse_graph->Draw("ALP");
        inverse_graph->GetXaxis()->SetRangeUser(get<0>(result) - 6, get<0>(result) + 3);
        fafter->SetRange(get<0>(result) - 6, get<0>(result) + 3);
        fafter->DrawCopy("SAME");
        TLatex *l1 = new TLatex();
        l1->DrawLatexNDC(0.2, 0.85, Form("V_{bd} = %.2f #pm %.2f", get<0>(result), get<1>(result)));
        //
        TString basedir = get_environment_variable("SIPM4EIC_WORKDIR");
        basedir += TString("/plots/vbdcheck/");
        gROOT->ProcessLine(Form(".! mkdir -p %s", basedir.Data()));
        c1->SaveAs(basedir + gTarget->GetName() + TString(".pdf"));
        delete c1;
    }
    delete log_graph;
    delete deriv_graph;
    delete inverse_graph;
    return result;
}
//
std::pair<float, float>
measure_breakdown_3(TGraphErrors *gTarget, float vbd_guess, bool graphics = false)
{
    std::pair<float, float> result = {-1, -1};
    auto log_graph = graphutils::log(gTarget);
    auto deriv_graph = graphutils::derivate(log_graph);
    auto deriv2_graph = graphutils::derivate(deriv_graph);
    TF1 *fpeak = new TF1("fpeak", "gaus(0)");
    deriv2_graph->Fit(fpeak, "Q", "SAME", vbd_guess - 0.7, vbd_guess + 0.7);
    get<0>(result) = fpeak->GetParameter(1);
    get<1>(result) = fpeak->GetParError(1);
    delete log_graph;
    delete deriv_graph;
    // delete deriv2_graph;
    return result;
}
//
std::pair<float, float>
measure_breakdown_4(TGraphErrors *gTarget, float vbd_guess, bool graphics = false)
{
    std::pair<float, float> result = {-1, -1};
    auto log_graph = graphutils::log(gTarget);
    //  Calculate a guess for the Vbd and have the measured baseline value
    auto vbd_guess_baseline = find_vbd_guess(log_graph);
    vbd_guess = get<0>(vbd_guess_baseline);
    auto mean_baseline = get<1>(vbd_guess_baseline);
    TF1 *fbefore = new TF1("fbefore", "[0]+x*[1]");
    TF1 *fafter = new TF1("fafter", "[0]+x*[1]+x*x*[2]");
    fbefore->SetParameter(0, mean_baseline);
    fbefore->SetParameter(1, 0.);
    fafter->SetParameter(0, 0.);
    fafter->SetParameter(1, 4);
    fafter->SetParLimits(2, -1.e3, 0.);
    fafter->SetParameter(2, -0.03);
    log_graph->Fit(fbefore, "Q", "", vbd_guess - 5, vbd_guess - 0.3);
    log_graph->Fit(fafter, "Q", "", vbd_guess + 0.3, vbd_guess + 2.3);
    auto intercept_before = fbefore->GetParameter(0);
    auto angularc_before = fbefore->GetParameter(1);
    auto c0_after = fafter->GetParameter(0);
    auto c1_after = fafter->GetParameter(1);
    auto c2_after = fafter->GetParameter(2);
    auto delta = sqrt((c1_after - angularc_before) * (c1_after - angularc_before) - 4 * c2_after * (c0_after - intercept_before));
    auto intercept_before_econtrib = fbefore->GetParError(0) / (delta);
    auto angularc_before_econtrib = fbefore->GetParError(1) * ((c1_after - angularc_before) / (delta) + 1) / (2 * c2_after);
    auto c0_after_econtrib = fafter->GetParError(0) / (delta);
    auto c1_after_econtrib = -fafter->GetParError(1) * ((c1_after - angularc_before) / (delta) + 1) / (2 * c2_after);
    auto c2_after_econtrib = fafter->GetParError(2) * ((c0_after - intercept_before) / (c2_after * delta) - (-delta - c1_after + angularc_before) / (2 * c2_after * c2_after));
    auto intercept_val = get<0>(result) = (-(c1_after - angularc_before) + delta) / (2 * c2_after);
    auto intercept_err = get<1>(result) = sqrt(intercept_before_econtrib * intercept_before_econtrib + angularc_before_econtrib * angularc_before_econtrib + c0_after_econtrib * c0_after_econtrib + c1_after_econtrib * c1_after_econtrib + c2_after_econtrib * c2_after_econtrib);
    TF1 *ffull = new TF1("ffull", "(0.5-0.5*TMath::Erf(((x-[0])/[1])))*([2]+x*[3])+(0.5+0.5*TMath::Erf(((x-[0])/[1])))*([4]+x*[5]+x*x*[6])");
    ffull->SetParameter(0, get<0>(result) + 0.2);
    ffull->SetParLimits(0, intercept_val - 0.3, intercept_val + 5.);
    ffull->SetParameter(1, 0.2);
    ffull->SetParLimits(1, 0.0001, 0.3);
    ffull->FixParameter(2, fbefore->GetParameter(0));
    ffull->FixParameter(3, fbefore->GetParameter(1));
    ffull->FixParameter(4, fafter->GetParameter(0));
    ffull->FixParameter(5, fafter->GetParameter(1));
    ffull->FixParameter(6, fafter->GetParameter(2));
    log_graph->Fit(ffull, "Q", "SAME", get<0>(result) - 4, get<0>(result) + 4.);
    // ffull->SetParLimits(0, intercept_val + ffull->GetParameter(1), intercept_val + 5.);
    // log_graph->Fit(ffull, "RQ", "SAME", get<0>(result) - 4, get<0>(result) + 3);
    get<0>(result) = ffull->GetParameter(0) - ffull->GetParameter(1);
    get<1>(result) = sqrt(ffull->GetParError(0) * ffull->GetParError(0) + ffull->GetParameter(1) * ffull->GetParameter(1));
    if (graphics)
    {
        TCanvas *c1 = new TCanvas(gTarget->GetName(), "c1", 600, 500);
        log_graph->Draw("ALP");
        log_graph->GetXaxis()->SetRangeUser(get<0>(result) - 6, get<0>(result) + 3);
        auto max = ffull->Eval(get<0>(result) + 2);
        max = max > 0 ? max * 1.1 : max * 0.9;
        log_graph->SetMaximum(max);
        fbefore->SetRange(get<0>(result) - 5, get<0>(result) + 5);
        fbefore->SetLineColor(kBlue);
        fbefore->DrawCopy("SAME");
        fafter->SetRange(get<0>(result) - 5, get<0>(result) + 5);
        fafter->SetLineColor(kGreen);
        fafter->DrawCopy("SAME");
        ffull->SetRange(get<0>(result) - 5, get<0>(result) + 5);
        ffull->DrawCopy("SAME");
        TLatex *l1 = new TLatex();
        l1->DrawLatexNDC(0.2, 0.85, Form("V_{bd} = %.2f #pm %.2f", get<0>(result), get<1>(result)));
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
std::pair<float, float>
measure_breakdown_5(TGraphErrors *gTarget, float vbd_guess, bool graphics = false)
{
    std::pair<float, float> result = {-1, -1};
    auto log_graph = graphutils::log(gTarget);
    TF1 *fbefore = new TF1("fbefore", "[0]+x*[1]");
    TF1 *fafter = new TF1("fafter", "[0]*(x-[1])^([2])");
    fafter->SetParameter(1, vbd_guess);
    log_graph->Fit(fbefore, "MERQ", "", vbd_guess - 5, vbd_guess - 0.2);
    log_graph->Fit(fafter, "MERQ", "", vbd_guess + 0.2, vbd_guess + 5);
    TF1 *ffull = new TF1("ffull", "([1]+x*[2])+(x>[0])*(+[3]*(x-[4])^([5]))");
    ffull->SetParameter(0, vbd_guess);
    ffull->FixParameter(1, fbefore->GetParameter(0));
    ffull->FixParameter(2, fbefore->GetParameter(1));
    ffull->FixParameter(3, fafter->GetParameter(0));
    ffull->FixParameter(4, fafter->GetParameter(1));
    ffull->FixParameter(5, fafter->GetParameter(2));
    log_graph->Fit(ffull, "MERQ", "SAME", vbd_guess - 3, vbd_guess + 3);
    get<0>(result) = 48; // ffull->GetParameter(0) - ffull->GetParameter(1);
    get<1>(result) = sqrt(ffull->GetParError(0) * ffull->GetParError(0) + ffull->GetParameter(1) * ffull->GetParameter(1));
    if (graphics)
    {
        TCanvas *c1 = new TCanvas();
        log_graph->Draw("ALP");
        log_graph->GetXaxis()->SetRangeUser(get<0>(result) - 2, get<0>(result) + 2);
        fbefore->SetRange(get<0>(result) - 5, get<0>(result) + 5);
        fbefore->SetLineColor(kBlue);
        fbefore->DrawCopy("SAME");
        fafter->SetRange(get<0>(result) - 5, get<0>(result) + 5);
        fafter->SetLineColor(kGreen);
        fafter->DrawCopy("SAME");
        ffull->SetRange(get<0>(result) - 5, get<0>(result) + 5);
        ffull->DrawCopy("SAME");
        TLatex *l1 = new TLatex();
        l1->DrawLatexNDC(0.2, 0.85, Form("V_{bd} = %.2f #pm %.2f", get<0>(result), get<1>(result)));
        //
        TString basedir = get_environment_variable("SIPM4EIC_WORKDIR");
        basedir += TString("/plots/vbdcheck/");
        gROOT->ProcessLine(Form(".! mkdir -p %s", basedir.Data()));
        c1->SaveAs(basedir + gTarget->GetName() + TString(".pdf"));
        delete c1;
    }
    delete log_graph;
    get<0>(result) = ffull->GetParameter(0);
    get<1>(result) = sqrt(ffull->GetParError(0) * ffull->GetParError(0));
    return result;
}
//
std::pair<float, float>
measure_breakdown_6(TGraphErrors *gTarget, float vbd_guess, bool graphics = false)
{
    std::pair<float, float> result = {-1, -1};
    auto current_graph = (TGraphErrors *)gTarget->Clone("tmp");
    TF1 *fprebefore = new TF1("fprebefore", "[0]");
    TF1 *fbefore = new TF1("fbefore", "[0]+x*[1]");
    TF1 *fafter = new TF1("fafter", "(TMath::Sign([3],(x-[0])))*(TMath::Power(TMath::Abs(x-[0])/([1]),[2]))");
    fbefore->FixParameter(0, fprebefore->GetParameter(0));
    fafter->FixParameter(0, vbd_guess);
    fafter->SetParLimits(1, 1.e-12, 1.e1);
    fafter->SetParameter(1, .7);
    fafter->SetParLimits(2, 1.e-12, 1.e1);
    fafter->SetParameter(2, .1);
    fafter->SetParameter(3, fprebefore->GetParameter(0));
    fafter->FixParameter(0, vbd_guess);
    fafter->SetParameter(1, 7.17717e-01);
    fafter->SetParameter(2, 1.);
    fafter->SetParameter(3, 4.59948e-09);
    current_graph->Fit(fbefore, "MERQ", "", vbd_guess - 5, vbd_guess - 0.2);
    current_graph->Fit(fafter, "MERQ", "", vbd_guess + 0.2, vbd_guess + 5);
    TF1 *ffull = new TF1("ffull", "(0.5-0.5*TMath::Erf((x-[0])/[1]))*([2]+x*[3])+(0.5+0.5*TMath::Erf(((x-[0])/[1])))*([2]+x*[3]+(TMath::Sign([7],(x-[4])))*(TMath::Power(TMath::Abs(x-[4])/([5]),[6])))");
    ffull->SetParameter(0, vbd_guess);
    ffull->SetParameter(1, 0.2);
    ffull->SetParLimits(1, 0.0001, 0.3);
    ffull->FixParameter(2, fbefore->GetParameter(0));
    ffull->FixParameter(3, fbefore->GetParameter(1));
    ffull->FixParameter(4, fafter->GetParameter(0));
    ffull->FixParameter(5, fafter->GetParameter(1));
    ffull->FixParameter(6, fafter->GetParameter(2));
    ffull->FixParameter(7, fafter->GetParameter(3));
    current_graph->Fit(ffull, "MERQ", "SAME", vbd_guess - 5, vbd_guess + 5);
    get<0>(result) = ffull->GetParameter(0) - ffull->GetParameter(1);
    get<1>(result) = sqrt(ffull->GetParError(0) * ffull->GetParError(0) + ffull->GetParameter(1) * ffull->GetParameter(1));
    if (graphics)
    {
        TCanvas *c1 = new TCanvas();
        gPad->SetLogy();
        current_graph->Draw("ALP");
        current_graph->GetXaxis()->SetRangeUser(get<0>(result) - 5, get<0>(result) + 5);
        fbefore->SetRange(get<0>(result) - 5, get<0>(result) + 5);
        fbefore->SetLineColor(kBlue);
        fbefore->DrawCopy("SAME");
        fafter->SetRange(get<0>(result) - 5, get<0>(result) + 5);
        fafter->SetLineColor(kGreen);
        fafter->DrawCopy("SAME");
        ffull->SetRange(get<0>(result) - 5, get<0>(result) + 5);
        ffull->DrawCopy("SAME");
        TLatex *l1 = new TLatex();
        l1->DrawLatexNDC(0.2, 0.85, Form("V_{bd} = %.2f #pm %.2f", get<0>(result), get<1>(result)));
        //
        TString basedir = get_environment_variable("SIPM4EIC_WORKDIR");
        basedir += TString("/plots/vbdcheck/");
        gROOT->ProcessLine(Form(".! mkdir -p %s", basedir.Data()));
        c1->SaveAs(basedir + gTarget->GetName() + TString(".pdf"));
        delete c1;
    }
    delete current_graph;
    get<0>(result) = ffull->GetParameter(0);
    get<1>(result) = sqrt(ffull->GetParError(0) * ffull->GetParError(0));
    return result;
}
//
std::pair<float, float>
measure_breakdown_100(TGraphErrors *gTarget, float vbd_guess, bool graphics = false)
{
    // Generate result pair
    std::pair<float, float> result = {-1, -1};
    // Geenrate local copy of graph
    auto current_graph = (TGraphErrors *)gTarget->Clone();
    // Linear fit to find with extrapolation the crossing of X-axis
    TF1 *fafter = new TF1("fafter", "(x-[0])/[1]");
    TF1 *ffull = new TF1("ffull", "(x-[0])/[1]+[2]/(x-[3])");
    ffull->SetParLimits(2, 0., 1.e8);
    // Function prepping
    // After fitting
    fafter->FixParameter(0, vbd_guess);
    fafter->SetParameter(1, 1. / (current_graph->Eval(vbd_guess + 3) - current_graph->Eval(vbd_guess + 2)));
    current_graph->Fit(fafter, "MEQ", "", vbd_guess + 2.0, vbd_guess + 2.5);
    fafter->SetParLimits(0, vbd_guess - 3, vbd_guess + 3);
    current_graph->Fit(fafter, "MEQ", "", vbd_guess + 2.0, vbd_guess + 2.5);
    // Full fitting
    ffull->FixParameter(0, fafter->GetParameter(0));
    ffull->FixParameter(1, fafter->GetParameter(1));
    ffull->SetParameter(2, 1.e+05);
    ffull->FixParameter(3, fafter->GetParameter(0) + 0.5);
    current_graph->Fit(ffull, "MEQ", "", vbd_guess + 0.5, vbd_guess + 2.5);
    ffull->SetParLimits(0, vbd_guess - 3, vbd_guess + 3);
    ffull->ReleaseParameter(1);
    ffull->SetParLimits(3, vbd_guess - 3, vbd_guess + 3);
    current_graph->Fit(ffull, "MEQ", "", vbd_guess + 0.5, vbd_guess + 2.5);
    // Set Parameters for show
    fafter->SetParameter(0, ffull->GetParameter(0));
    fafter->SetParameter(1, ffull->GetParameter(1));
    // Save result
    get<0>(result) = ffull->GetParameter(0);
    get<1>(result) = ffull->GetParError(0);
    // Plot if requested
    if (graphics)
    {
        TCanvas *c1 = new TCanvas();
        gPad->SetLogy();
        current_graph->Draw("ALP");
        current_graph->GetXaxis()->SetRangeUser(get<0>(result) - 5, get<0>(result) + 5);
        fafter->SetRange(get<0>(result) - 5, get<0>(result) + 5);
        fafter->SetLineColor(kBlue);
        fafter->DrawCopy("SAME");
        ffull->SetRange(get<0>(result) - 5, get<0>(result) + 5);
        ffull->SetLineColor(kRed);
        ffull->DrawCopy("SAME");
        // Declare found Vbd value
        TLatex *l1 = new TLatex();
        l1->DrawLatexNDC(0.2, 0.85, Form("V_{bd} = %.2f #pm %.2f", get<0>(result), get<1>(result)));
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
std::pair<float, float>
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
std::pair<float, float>
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
std::pair<float, float>
overvoltage_ratio(TGraphErrors *gNumerator, TGraphErrors *gDenominator, float guess_vbd_num, float guess_vbd_den, float _xcoordinate)
{
    std::pair<float, float> result;
    std::pair<float, float> numerator;
    std::pair<float, float> denominator;
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
    get<0>(result) = get<0>(numerator) / get<0>(denominator);
    get<1>(result) = get<0>(result) * ((get<1>(numerator) / get<0>(numerator)) * (get<1>(numerator) / get<0>(numerator)) + (get<1>(denominator) / get<0>(denominator)) * (get<1>(denominator) / get<0>(denominator)));
    return result;
}
//
template <bool recalulate_vbd = true>
TGraphErrors *
overvoltage_ratio(TGraphErrors *gNumerator, TGraphErrors *gDenominator, float guess_vbd_num, float guess_vbd_den, std::vector<float> _xcoordinates)
{
    TGraphErrors *result = new TGraphErrors();
    std::pair<float, float> current_ratio;
    std::pair<float, float> numerator;
    std::pair<float, float> denominator;
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
    std::pair<float, float> current_ratio;
    std::pair<float, float> numerator;
    std::pair<float, float> denominator;
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
}