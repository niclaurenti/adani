#include "adani/Exceptions.h"

//==========================================================================================//
//  Definition of colors
//------------------------------------------------------------------------------------------//

#define RESET "\033[0m"
#define RED "\033[31m"
#define GREEN "\033[32m"
#define YELLOW "\033[33m"
#define BLUE "\033[34m"

//==========================================================================================//
//  UnknownException: constructor
//------------------------------------------------------------------------------------------//

UnknownException::UnknownException(
    const string &reason, const string &function, const int &line
) {
    message_ = "UnknownException in " + function + " at line " + to_string(line) + ": " + reason;
}

//==========================================================================================//
//  UnknownException: redefinition of what method
//------------------------------------------------------------------------------------------//

const char *UnknownException::what() const noexcept {
    return message_.c_str();
}

//==========================================================================================//
//  UnknownException: runtime_error
//------------------------------------------------------------------------------------------//

void UnknownException::runtime_error() const {
    cerr << RED << "Error! " << RESET;
    cerr << (*this).what() << endl;
    cerr << "Exiting program!" << endl;
    exit(-1);
}

//==========================================================================================//
//  UnknownException: warning
//------------------------------------------------------------------------------------------//

void UnknownException::warning() const {
    cerr << YELLOW << "Warning! " << RESET;
    cerr << (*this).what() << endl;
}

//==========================================================================================//
//  NotImplementedException: constructor
//------------------------------------------------------------------------------------------//

NotImplementedException::NotImplementedException(
    const string &reason, const string &function, const int &line
) {
    message_ = "NotImplementedException in " + function + " at line " + to_string(line) + ": " + reason;
}

//==========================================================================================//
//  NotKnownException: constructor
//------------------------------------------------------------------------------------------//

NotKnownException::NotKnownException(
    const string &reason, const string &function, const int &line
) {
    message_ = "NotKnownException in " + function + " at line " + to_string(line) + ": " + reason;
}

//==========================================================================================//
//  NotPresentException: constructor
//------------------------------------------------------------------------------------------//

NotPresentException::NotPresentException(
    const string &reason, const string &function, const int &line
) {
    message_ = "NotPresentException in " + function + " at line " + to_string(line) + ": " + reason;
}

//==========================================================================================//
//  NotValidException: constructor
//------------------------------------------------------------------------------------------//

NotValidException::NotValidException(
    const string &reason, const string &function, const int &line
) {
    message_ = "NotValidException in " + function + " at line " + to_string(line) + ": " + reason;
}

//==========================================================================================//
//  UnexpectedException: constructor
//------------------------------------------------------------------------------------------//

UnexpectedException::UnexpectedException(
    const string &reason, const string &function, const int &line
) {
    message_ = "UnexpectedException in " + function + " at line " + to_string(line) + ": " + reason;
}
