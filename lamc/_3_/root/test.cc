#include <set>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <filesystem>

namespace fs = std::filesystem;

std::vector<std::string> split (const std::string &s, char delim) {
    std::vector<std::string> result;
    std::stringstream ss (s);
    std::string item;

    while (getline (ss, item, delim)) {
        result.push_back (item);
    }

    return result;
}

int main() {



  std::string path_name = fs::current_path();

  //--- filenames are unique so we can use a set
  std::set<fs::path> sorted_by_name;

  for (auto &entry : fs::directory_iterator(path_name))
    sorted_by_name.insert(entry.path());

  //--- print the files sorted by filename
  for (auto &filename : sorted_by_name)
    std::cout << filename.c_str() << std::endl;


  return 0;
}
