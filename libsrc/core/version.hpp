
namespace ngcore
{
  class VersionInfo
  {
  private:
    size_t mayor, minor, date, commit_offset;
    std::string git_hash;
  public:
    VersionInfo() : mayor(0), minor(0), date(0), commit_offset(0), git_hash("") {}
    VersionInfo(std::string vstring)
    {
      minor = date = commit_offset = 0;
      git_hash = "";
      if(vstring.substr(0,1) == "v")
        vstring = vstring.substr(1,vstring.size()-1);
      auto dot = vstring.find(".");
      mayor = std::stoi(vstring.substr(0,dot));
      if(dot == size_t(-1)) vstring = "";
      else vstring = vstring.substr(dot+1, vstring.size()-dot-1);
      if(vstring.size())
        {
          dot = vstring.find(".");
          minor = std::stoi(vstring.substr(0,dot));
          if (dot == size_t(-1)) vstring = "";
          else vstring = vstring.substr(dot+1, vstring.size()-dot-1);
          if(vstring.size())
            {
              dot = vstring.find("-");
              date = std::stoi(vstring.substr(0,dot));
              if(dot == size_t(-1)) vstring = "";
              else vstring = vstring.substr(dot+1,vstring.size()-dot-1);
              if(vstring.size())
                {
                  dot = vstring.find("-");
                  commit_offset = std::stoi(vstring.substr(0,dot));
                  if(dot == size_t(-1)) vstring = "";
                  else vstring = vstring.substr(dot+1, vstring.size()-dot-1);
                  if(vstring.size())
                    git_hash = vstring;
                }
            }
        }
    }
    VersionInfo(const char* cstr) : VersionInfo(std::string(cstr)) { }

    std::string to_string() const
    { std::string vstring = "v" + std::to_string(mayor);
      if(minor || date || commit_offset || git_hash.size())
        {
          vstring += "." + std::to_string(minor);
          if(date || commit_offset || git_hash.size())
            {
              vstring += "." + std::to_string(date);
              if(commit_offset || git_hash.size())
                {
                  vstring += "-" + std::to_string(commit_offset);
                  if(git_hash.size())
                    vstring += "-" + git_hash;
                }
            }
        }
      return vstring;
    }
    bool operator <(const VersionInfo& other) const
    {
      return std::tie(mayor, minor, date, commit_offset) <
        std::tie(other.mayor, other.minor, other.date, other.commit_offset);
    }
    bool operator ==(const VersionInfo& other) const
    {
      return mayor == other.mayor && minor == other.minor && date == other.date
        && commit_offset == other.commit_offset;
    }
    bool operator >(const VersionInfo& other) const { return other < (*this); }
    bool operator <=(const VersionInfo& other) const { return !((*this) > other); }
    bool operator >=(const VersionInfo& other) const { return !((*this) < other); }

    void DoArchive(Archive& ar)
    {
      ar & mayor & minor & date & commit_offset & git_hash;
    }
  };
}
