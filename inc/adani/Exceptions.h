/*
 * =====================================================================================
 *
 *       Filename:  Exceptions.h
 *
 *         Author:  Daniele Adani
 *
 *    Description:  Noi parliamo di idee che danno animo a un gioco che poi è
 *                  quello che amiamo, ovvero il calcio
 *
 *  In this file there are the custom exceptions
 *
 * =====================================================================================
 */
#ifndef Exceptions_h
#define Exceptions_h

#include <iostream>
#include <stdexcept>

using namespace std;

//==========================================================================================//
//  class UnknownException: base class for all the exceptions
//------------------------------------------------------------------------------------------//

class UnknownException : public exception {
    public:
        UnknownException(){};
        UnknownException(const string &reason, const string &function, const int &line);

        const char *what() const noexcept override;
        void runtime_error() const;
        void warning() const;

    protected:
        string message_;
};

//==========================================================================================//
//  class NotImplementedException: when something is not implemented
//------------------------------------------------------------------------------------------//

class NotImplementedException : public UnknownException {
    public:
        NotImplementedException(const string &reason, const string &function, const int &line);
};

//==========================================================================================//
//  class NotKnownException: when something is still not known
//------------------------------------------------------------------------------------------//

class NotKnownException : public UnknownException {
    public:
        NotKnownException(const string &reason, const string &function, const int &line);
};

//==========================================================================================//
//  class NotPresentException: when something is not present, i.e. it is zero
//------------------------------------------------------------------------------------------//

class NotPresentException : public UnknownException {
    public:
        NotPresentException(const string &reason, const string &function, const int &line);
};

//==========================================================================================//
//  class NotValidException: passing invalid parameters
//------------------------------------------------------------------------------------------//

class NotValidException : public UnknownException {
    public:
        NotValidException(const string &reason, const string &function, const int &line);
};

//==========================================================================================//
//  class UnexpectedException: passing invalid parameters
//------------------------------------------------------------------------------------------//

class UnexpectedException : public UnknownException {
    public:
        UnexpectedException(const string &reason, const string &function, const int &line);
};

#endif
