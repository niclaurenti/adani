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
//  class AbstractException: base class for all the exceptions
//------------------------------------------------------------------------------------------//

class AbstractException : public exception {
    public:
        AbstractException(){};
        AbstractException(const string &reason, const string &function);

        const char *what() const noexcept override;
        void runtime_error() const;
        void warning() const;

    protected:
        string message_;
};

//==========================================================================================//
//  class NotImplementedException: when something is not implemented
//------------------------------------------------------------------------------------------//

class NotImplementedException : public AbstractException {
    public:
        NotImplementedException(const string &reason, const string &function);
};

//==========================================================================================//
//  class NotKnownException: when something is still not known
//------------------------------------------------------------------------------------------//

class NotKnownException : public AbstractException {
    public:
        NotKnownException(const string &reason, const string &function);
};

//==========================================================================================//
//  class NotPresentException: when something is not present, i.e. it is zero
//------------------------------------------------------------------------------------------//

class NotPresentException : public AbstractException {
    public:
        NotPresentException(const string &reason, const string &function);
};

//==========================================================================================//
//  class NotValidException: passing invalid parameters
//------------------------------------------------------------------------------------------//

class NotValidException : public AbstractException {
    public:
        NotValidException(const string &reason, const string &function);
};

//==========================================================================================//
//  class UnexpectedException: passing invalid parameters
//------------------------------------------------------------------------------------------//

class UnexpectedException : public AbstractException {
    public:
        UnexpectedException(const string &reason, const string &function);
};

#endif
