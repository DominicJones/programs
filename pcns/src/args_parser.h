// -*- C++ -*-
#pragma once

#include <string>
#include <vector>
#include <algorithm>

class args_parser_t
{
public:
  args_parser_t(int argc, char** argv)
  {
    for (auto it = 0; it != argc; ++it)
      _tokens.push_back(std::string(argv[it]));
  }

  auto size() const { return _tokens.size(); };
  auto arg(int const &it) const { return _tokens[it]; };

  auto get_option(const std::string &option) const
  {
    auto itr = std::find(_tokens.begin(), _tokens.end(), option);
    if (itr != _tokens.end() && ++itr != _tokens.end()) return std::string{*itr};
    return std::string{};
  }

  bool has_option(const std::string &option) const
  {
    return std::find(_tokens.begin(), _tokens.end(), option) != _tokens.end();
  }

private:
  std::vector<std::string> _tokens;
};
