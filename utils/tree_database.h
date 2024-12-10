#pragma once

#include "general_utility.h"

class tree_database
{
  std::unordered_map<int, std::string> name_database;                  // Mapping from node ID to display name
  std::unordered_map<std::string, std::vector<int>> omonymes_database; // Mapping from display names to node IDs
  std::map<int, std::vector<int>> database;                            // Tree structure with parent ID -> vector of child IDs
  int current_id = 0;                                                  // To generate unique IDs

public:
  //  Declaration
  //  --- I/O
  int add_node(int parent_id, const std::string &child_name);
  int inline add_node(const std::string &child_name) { return add_node(-1, child_name); };
  template <typename... Args>
  int add_node(const Args &...ancestors, const std::string &child_name);
  int find_node(int parent_id, const std::string &name);
  int find_node(int parent_id, int depth, const std::string &name);
  int find_or_create_node(int parent_id, const std::string &name);
  template <typename... Args>
  int find_or_create_node(int parent_id, const Args &...ancestors);
  //  --- Getters
  std::vector<int> get_children_ids(int parent_id);
  std::vector<std::pair<int, std::string>> get_children(int parent_id) { return get_ids_name(get_children_ids(parent_id)); }
  std::vector<int> get_children_ids(int parent_id, int depth);
  std::vector<std::pair<int, std::string>> get_children(int parent_id, int depth) { return get_ids_name(get_children_ids(parent_id, depth)); }
  std::vector<int> get_children_ids(int parent_id, int depth, std::vector<std::pair<int, std::string>> filters);
  std::vector<std::pair<int, std::string>> get_children(int parent_id, int depth, std::vector<std::pair<int, std::string>> filters) { return get_ids_name(get_children_ids(parent_id, depth, filters)); }
  std::vector<int> get_name_ids(std::string name) { return omonymes_database[name]; }
  std::pair<int, std::string> get_id_name(int id) { return {id, name_database[id]}; }
  std::vector<std::pair<int, std::string>> get_ids_name(std::vector<int> id_list);
  int get_ancestor_id(int child_id, int depth);
  std::pair<int, std::string> get_ancestor(int child_id, int depth) { return get_id_name(get_ancestor_id(child_id, depth)); }
  //  --- Helpers
  bool is_name_unique(std::string name) { return get_name_ids(name).size() == 1; }
  // --- Print function
  void print_tree(int node_id, int depth);
  void inline print_tree() { return print_tree(-1, -1); }
};

//  Implementation
//  --- I/O
int tree_database::add_node(int parent_id, const std::string &child_name)
{
  int child_id = current_id++; // Generate unique ID for the child
  name_database[child_id] = child_name;
  omonymes_database[child_name].push_back(child_id);
  database[parent_id].push_back(child_id);
  return child_id;
}
template <typename... Args>
int tree_database::add_node(const Args &...ancestors, const std::string &child_name)
{
  // Check all ancestors are string
  static_assert((std::is_same_v<Args, std::string> && ...), "[ERROR][database::tree_database::add_node] All ancestors must be of type std::string");

  std::vector<std::string> nodes = {ancestors...};
  int parent_id = -1;

  // Traverse through ancestors, ensuring each exists or is created
  for (const auto &node : nodes)
    parent_id = find_or_create_node(parent_id, node);

  // Now add the final child node under the last parent
  return find_or_create_node(parent_id, child_name);
}
int tree_database::find_node(int parent_id, const std::string &name)
{
  auto node_children = get_children(parent_id);
  for (auto id : node_children)
    if (id.second == name)
      return id.first; // Node exists, return its ID
  return -1;
}
int tree_database::find_node(int parent_id, int depth, const std::string &name)
{
  auto node_children = get_children(parent_id, depth);
  for (auto id : node_children)
    if (name_database[id.first] == name)
      return id.first; // Node exists, return its ID
  return -1;
}
int tree_database::find_or_create_node(int parent_id, const std::string &name)
{
  auto found_node = find_node(parent_id, name);
  if (found_node >= 0)
    return found_node;

  // Node doesn't exist, create a new one
  return add_node(parent_id, name);
}
template <typename... Args>
int tree_database::find_or_create_node(int parent_id, const Args &...ancestors)
{
  static_assert((std::is_same_v<Args, std::string> && ...), "[ERROR][tree_database::find_or_create_node] All ancestors must be of type std::string");

  std::vector<std::string> nodes = {ancestors...};

  // Traverse through ancestors, ensuring each exists or is created
  for (const auto &node : nodes)
  { // cout << node << endl;
    parent_id = find_or_create_node(parent_id, node);
  }

  return parent_id; // Return the last node ID created or found
}
//  --- Getters
std::vector<int> tree_database::get_children_ids(int parent_id)
{
  std::vector<int> children;
  for (auto current_child : database[parent_id])
    children.push_back(current_child);
  return children;
}
std::vector<int> tree_database::get_children_ids(int parent_id, int depth)
{
  std::vector<int> result;
  std::vector<int> prev_gen;
  if (depth <= 0)
    return result;
  auto first_degree_children = get_children_ids(parent_id);
  if (depth == 1)
    return first_degree_children;
  prev_gen = first_degree_children;
  for (auto i_depth = 1; i_depth < depth; i_depth++)
  {
    result.clear();
    for (auto current_parent : prev_gen)
      result = utility::merge(result, get_children_ids(current_parent));
    prev_gen = result;
  }
  return result;
}
std::vector<int> tree_database::get_children_ids(int parent_id, int depth, std::vector<std::pair<int, std::string>> filters)
{
  std::vector<std::pair<int, std::string>> result;

  //  Trivial case, no reasonable filter is applicable
  if (depth < 2)
    return get_children_ids(parent_id, depth);

  // Filter out any filters with depths below 1 or greater than the specified depth
  filters.erase(std::remove_if(filters.begin(), filters.end(), [depth](const std::pair<int, std::string> &filter)
                               { return filter.first < 1 || filter.first > depth; }),
                filters.end());

  if (filters.size() == 0)
    return get_children_ids(parent_id, depth);

  // Map to group filters by depth
  std::unordered_map<int, std::vector<std::string>> filters_names_per_depth;
  std::unordered_map<int, std::vector<int>> filters_ids_per_depth;
  for (const auto &filter : filters)
  {
    int filter_depth = filter.first;
    std::string filter_name = filter.second;
    filters_names_per_depth[filter_depth].push_back(filter_name);
  }
  std::vector<int> current_parent_ids = {parent_id};
  std::vector<int> current_children_ids = {};
  for (auto i_depth = 1; i_depth <= depth; i_depth++)
  {
    for (auto current_parent : current_parent_ids)
      if (filters_names_per_depth[i_depth].size() > 0)
        for (auto current_filter_name : filters_names_per_depth[i_depth])
        {
          auto found_node = find_node(current_parent, current_filter_name);
          if (found_node >= 0)
            current_children_ids.push_back(found_node);
        }
      else
        current_children_ids = utility::merge(current_children_ids, get_children_ids(current_parent));
    current_parent_ids = current_children_ids;
    current_children_ids.clear();
  }
  return current_parent_ids;
}
std::vector<std::pair<int, std::string>> tree_database::get_ids_name(std::vector<int> id_list)
{
  std::vector<std::pair<int, std::string>> result;
  for (auto current_id : id_list)
    result.push_back({current_id, name_database[current_id]});
  return result;
}
int tree_database::get_ancestor_id(int child_id, int depth)
{
  if (depth == 0)
    return child_id;
  std::pair<int, std::string> result;
  for (auto [parent_id, children_ids] : database)
  {
    auto current_child = std::find(children_ids.begin(), children_ids.end(), child_id);
    if (current_child != children_ids.end())
      return get_ancestor_id(parent_id, depth - 1);
  }
  return -1;
}
// --- Print function
void tree_database::print_tree(int node_id, int depth)
{
  // Indent based on the depth of the node in the tree
  std::cout << "|- ";
  if (node_id == depth && depth == -1)
    std::cout << "Printing database";
  for (int i = 0; i < depth; ++i)
    std::cout << " - "; // Characters per level of depth

  // Print the current node's name
  std::cout << name_database[node_id] << std::endl;

  // If the node has children, recursively print them
  if (database.find(node_id) != database.end())
    for (int child_id : database[node_id])
      print_tree(child_id, depth + 1); // Recursive call for each child
}
