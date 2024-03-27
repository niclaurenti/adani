/*
 * =====================================================================================
 *
 *       Filename:  SplittingFunction.h
 *
 *    Description:  Header file for the
 * SplittingFunction.cc file.
 *
 *         Author:  A livello di servilismo come siamo
 * messi?
 *
 *  In this file there are the splitting functions.
 *
 * =====================================================================================
 */

#ifndef Split
#define Split

//==========================================================================================//
//  class AbstractSplittingFunction
//------------------------------------------------------------------------------------------//

class AbstractSplittingFunction {
    public:
        AbstractSplittingFunction() { mult_factor_ = 1.; };
        virtual ~AbstractSplittingFunction() = 0;

        // Components of the Splitting Function
        virtual double Regular(double x, int nf) const = 0;
        virtual double Singular(double x, int nf) const = 0;
        virtual double Local(int nf) const = 0;
        // Integral from 0 to x of the Singular part
        virtual double SingularIntegrated(double x, int nf) const = 0;

        double GetMultFact() const { return mult_factor_; };
        void SetMultFact(const double &mult_factor) {
            mult_factor_ = mult_factor;
        };

    private:
        double mult_factor_;
};

//==========================================================================================//
//  class SplittingFunction
//------------------------------------------------------------------------------------------//

class SplittingFunction : public AbstractSplittingFunction {
    public:
        SplittingFunction(
            const int &order, const char &entry1, const char &entry2
        );
        ~SplittingFunction() override{};

        // overloading operators
        SplittingFunction operator*(const double &rhs) const;
        friend SplittingFunction
        operator*(const double &lhs, const SplittingFunction &rhs);
        SplittingFunction operator/(const double &rhs) const;

        // Components of the Splitting Function
        double Regular(double x, int nf) const override;
        double Singular(double x, int nf) const override;
        double Local(int nf) const override;
        // Integral from 0 to x of the Singular part
        double SingularIntegrated(double x, int nf) const override;

        void SetFunctions();

        // get methods
        double GetOrder() const { return order_; };
        char GetEntry1() const { return entry1_; };
        char GetEntry2() const { return entry2_; };

    private:
        int order_;
        char entry1_;
        char entry2_;

        double (SplittingFunction::*reg_)(double, int) const;
        double (SplittingFunction::*sing_)(double, int) const;
        double (SplittingFunction::*sing_int_)(double, int) const;
        double (SplittingFunction::*loc_)(int) const;

        //==========================================================================================//
        //                      Splitting functions O(as)
        //                      without color factors
        //------------------------------------------------------------------------------------------//

        double pgq(double x) const;
        double pqg(double x) const;
        double pggreg(double x) const;
        double pggsing(double x) const;

        //==========================================================================================//
        //                      Splitting functions O(as)
        //------------------------------------------------------------------------------------------//

        double Pgq0(double x, int /*nf*/) const;
        double Pqg0(double x, int nf) const;

        double Pgg0reg(double x, int /*nf*/) const;
        double Pgg0loc(int nf) const;
        double Pgg0sing(double x, int /*nf*/) const;
        double Pgg0sing_integrated(double x, int /*nf*/) const;

        double Pqq0reg(double x, int /*nf*/) const;
        double Pqq0loc(int /*nf*/) const;
        double Pqq0sing(double x, int /*nf*/) const;
        double Pqq0sing_integrated(double x, int /*nf*/) const;

        //==========================================================================================//
        //                      Splitting functions O(as^2)
        //------------------------------------------------------------------------------------------//

        double Pgq1(double x, int nf) const;

        double Pgg1reg(double x, int nf) const;
        double Pgg1sing(double x, int nf) const;
        double Pgg1sing_integrated(double x, int nf) const;
        double Pgg1loc(int nf) const;

        //==========================================================================================//
        //  Function needed to return a zero function
        //------------------------------------------------------------------------------------------//

        double ZeroFunction_x_nf(double /*x*/, int /*nf*/) const { return 0.; };
        double ZeroFunction_nf(int /*nf*/) const { return 0.; };
};

//==========================================================================================//
//  class ConvolutedSplittingFunctions: analytical convolution between two
//  Splitting Functions
//------------------------------------------------------------------------------------------//

class ConvolutedSplittingFunctions : public AbstractSplittingFunction {
    public:
        ConvolutedSplittingFunctions(
            const int &order1, const char &entry1, const char &entry2,
            const int &order2, const char &entry3, const char &entry4
        );
        ~ConvolutedSplittingFunctions() override{};

        // overloading operators
        ConvolutedSplittingFunctions operator*(const double &rhs) const;
        friend ConvolutedSplittingFunctions
        operator*(const double &lhs, const ConvolutedSplittingFunctions &rhs);
        ConvolutedSplittingFunctions operator/(const double &rhs) const;

        // Components of the Convoluted Splitting Function
        double Regular(double x, int nf) const override;
        double Singular(double /*x*/, int /*nf*/) const override { return 0.; };
        double Local(int /*nf*/) const override { return 0.; };
        // Integral from 0 to x of the Singular part
        double SingularIntegrated(double /*x*/, int /*nf*/) const override {
            return 0.;
        };

        void SetFunctions();

        // get methods
        double GetOrder1() const { return order1_; };
        char GetEntry1() const { return entry1_; };
        char GetEntry2() const { return entry2_; };
        double GetOrder2() const { return order2_; };
        char GetEntry3() const { return entry3_; };
        char GetEntry4() const { return entry4_; };

    private:
        int order1_;
        char entry1_;
        char entry2_;
        int order2_;
        char entry3_;
        char entry4_;

        double (ConvolutedSplittingFunctions::*reg_)(double, int) const;

        //==========================================================================================//
        //  Convolution between the splitting functions Pgg0/Pqq0
        //  and Pgq0
        //------------------------------------------------------------------------------------------//

        double Pgg0_x_Pgq0(double x, int nf) const;
        double Pqq0_x_Pgq0(double x, int /*nf*/) const;

        //==========================================================================================//
        //  Convolution between the splitting functions Pgq0 and Pqg0
        //------------------------------------------------------------------------------------------//

        double Pgq0_x_Pqg0(double x, int nf) const;

        //==========================================================================================//
        //  Function needed to return a zero function
        //------------------------------------------------------------------------------------------//

        double ZeroFunction_x_nf(double /*x*/, int /*nf*/) const { return 0.; };
        double ZeroFunction_nf(int /*nf*/) const { return 0.; };
};

//==========================================================================================//
//  class Delta
//------------------------------------------------------------------------------------------//

class Delta : public AbstractSplittingFunction {
    public:
        Delta() : AbstractSplittingFunction(){};
        ~Delta() override{};

        double Regular(double /*x*/, int /*nf*/) const override { return 0.; };
        double Singular(double /*x*/, int /*nf*/) const override { return 0.; };
        double Local(int /*nf*/) const override { return GetMultFact() * 1.; };
        double SingularIntegrated(double /*x*/, int /*nf*/) const override {
            return 0.;
        };

        Delta operator*(const double &rhs) const;
        friend Delta operator*(const double &lhs, const Delta &rhs);

        Delta operator/(const double &rhs) const;
};

#endif
