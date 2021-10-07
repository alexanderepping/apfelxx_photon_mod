//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/messages.h"
#include "apfel/version.h"

namespace apfel
{
  /**
   * @brief Verbosity enumerator
   *
   * LOW    = No informative and warning messages are displayed. Error messages are printed anyway.
   * MEDIUM = No informative messages are displayed. Warning and error messages are printed.
   * HIGH   = All messages are displayed.
   */
  enum verbosity: int {LOW, MEDIUM, HIGH};

  /**
   * @brief Current Verbosity level
   */
  static int VerbosityLevel = HIGH;

  /**
   * @brief Colour codes
   */
  enum ColourCode {red = 31, green = 32, yellow = 33, blue = 34, normal = 39};

  //_________________________________________________________________________
  void SetVerbosityLevel(int const& vl)
  {
    VerbosityLevel = vl;
  }

  //_________________________________________________________________________
  int GetVerbosityLevel()
  {
    return VerbosityLevel;
  }

  //_________________________________________________________________________
  void report(std::string const& what)
  {
    if (VerbosityLevel > MEDIUM)
      std::cout << what;
  }

  //_________________________________________________________________________
  void info(std::string const& tag, std::string const& what)
  {
    if (VerbosityLevel > MEDIUM)
      std::cout << "\033[" << ColourCode::blue << "m[apfel::" << tag << "] Info: " << what << "\033[" << ColourCode::normal << "m\n";
  }

  //_________________________________________________________________________
  void warning(std::string const& tag, std::string const& what)
  {
    if (VerbosityLevel > LOW)
      std::cout << "\033[" << ColourCode::yellow << "m[apfel::" << tag << "] Warning: " << what << "\033[" << ColourCode::normal << "m\n";
  }

  //_________________________________________________________________________
  std::string error(std::string const& tag, std::string const& what)
  {
    std::stringstream ss;
    ss << "\n\n\033[" << ColourCode::red << "m[apfel::" << tag << "] Error: " << what << "\033[" << ColourCode::normal << "m\n";
    return ss.str();
  }

  //_________________________________________________________________________
  void Banner()
  {
    if (VerbosityLevel > MEDIUM)
      {
        std::cout << "\033[1;" << ColourCode::red << "m\n";
        std::cout << "       _/_/_/   _/_/_/_/  _/_/_/_/  _/_/_/_/  _/\n";
        std::cout << "     _/    _/  _/    _/  _/        _/        _/        _/     _/\n";
        std::cout << "    _/_/_/_/  _/_/_/_/  _/_/_/    _/_/_/    _/      _/_/_/ _/_/_/\n";
        std::cout << "   _/    _/  _/        _/        _/        _/        _/     _/\n";
        std::cout << "  _/    _/  _/        _/        _/_/_/_/  _/_/_/_/\n";
        std::cout << "\033[0;" << ColourCode::green << "m\n";
        std::cout << "_____v" << VERSION << ": A PDF evolution library in C++, arXiv:1708.00911\n";
        std::cout << "     Author: V. Bertone\n";
        std::cout << "\033[" << ColourCode::normal << "m\n";
      }
  }
}
