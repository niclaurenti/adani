#include "adani/Exceptions.h"

//==========================================================================================//
//  Definition of colors
//------------------------------------------------------------------------------------------//

#define RESET   "\033[0m"
#define RED     "\033[31m"
#define GREEN   "\033[32m"
#define YELLOW  "\033[33m"
#define BLUE    "\033[34m"

//==========================================================================================//
//  AbstractException: constructor
//------------------------------------------------------------------------------------------//

AbstractException::AbstractException(const string& reason, const string& function) {
    message_ = "AbstractException in " + function + ": " + reason;
}

//==========================================================================================//
//  AbstractException: redefinition of what method
//------------------------------------------------------------------------------------------//

const char* AbstractException::what() const noexcept {
    return message_.c_str();
}

//==========================================================================================//
//  AbstractException: runtime_error
//------------------------------------------------------------------------------------------//

void AbstractException::runtime_error() const {
    cerr << RED << "Error! " << RESET;
    cerr << (*this).what() << endl;
    cerr << "Exiting program!" << endl;
    exit(-1);
}

//==========================================================================================//
//  AbstractException: warning
//------------------------------------------------------------------------------------------//

void AbstractException::warning() const {
    cerr << YELLOW << "Warning! " << RESET;
    cerr << (*this).what() << endl;
}

//==========================================================================================//
//  NotImplementedException: constructor
//------------------------------------------------------------------------------------------//

NotImplementedException::NotImplementedException(const string& reason, const string& function) {
    message_ = "NotImplementedException in " + function + ": " + reason;
}

//==========================================================================================//
//  NotKnownException: constructor
//------------------------------------------------------------------------------------------//

NotKnownException::NotKnownException(const string& reason, const string& function) {
    message_ = "NotKnownException in " + function + ": " + reason;
}

//==========================================================================================//
//  NotPresentException: constructor
//------------------------------------------------------------------------------------------//

NotPresentException::NotPresentException(const string& reason, const string& function) {
    message_ = "NotPresentException in " + function + ": " + reason;
}

//==========================================================================================//
//  NotValidException: constructor
//------------------------------------------------------------------------------------------//

NotValidException::NotValidException(const string& reason, const string& function) {
    message_ = "NotValidException in " + function + ": " + reason;
}
