/*
 * =====================================================================================
 *
 *       Filename:  SplittingFunctions.h
 *
 *    Description:  Header file for the
 * SplittingFunctions.cc file.
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

class Function {
    public:
        virtual double Regular(const double x, const int nf) const = 0;
        virtual double Singular(const double x, const int nf) const = 0;
        virtual double Local(const int nf) const = 0;
        virtual double SingularIntegrated(const double x, const int nf) const = 0;

    double GetMultFact() const {return mult_factor_;};

    private:
        double mult_factor_;
};

class SplittingFunction : public Function{
    public:
        SplittingFunction(const int& order, const char& entry1, const char& entry2) ;
        ~SplittingFunction() {} ;

        double Regular(const double x, const int nf) const ;
        double Singular(const double x, const int nf) const ;
        double Local(const int nf) const ;
        double SingularIntegrated(const double x, const int nf) const;

        // SplittingFunction operator+(const SplittingFunction& splitfunc) const;

        // get methods
        double GetOrder() const {return order_ ;} ;
        char GetEntry1() const {return entry1_;} ;
        char GetEntry2() const {return entry2_;} ;

    private:
        int order_;
        char entry1_;
        char entry2_;

        //==========================================================================================//
        //                      Splitting functions O(alpha_s)
        //                      without color factors
        //------------------------------------------------------------------------------------------//

        double pgq(const double x) const ;
        double pqg(const double x) const ;
        double pggreg(const double x) const ;
        double pggsing(const double x) const ;

        //==========================================================================================//
        //                      Splitting functions O(alpha_s)
        //------------------------------------------------------------------------------------------//

        double Pgq0(const double x) const ;
        double Pqg0(const double x, const int nf) const ;

        double Pgg0reg(const double x) const ;
        double Pgg0loc(const int nf) const ;
        double Pgg0sing(const double x) const ;
        double Pgg0sing_integrated(const double x) const;

        double Pqq0reg(const double x) const ;
        double Pqq0loc() const ;
        double Pqq0sing(const double x) const ;
        double Pqq0sing_integrated(const double x) const;

        //==========================================================================================//
        //                      Splitting functions O(alpha_s^2)
        //------------------------------------------------------------------------------------------//

        double Pgq1(const double x, const int nf) const ;

        double Pgg1reg(const double x, const int nf) const ;
        double Pgg1sing(const double x, const int nf) const ;
        double Pgg1sing_integrated(const double x, const int nf) const ;
        double Pgg1loc(const int nf) const ;

};

class ConvolutedSplittingFunctions : public SplittingFunction {
    public:
        ConvolutedSplittingFunctions(const int& order, const char& entry1, const char& entry2, const char& entry3);

        char GetEntry3() const {return entry3_;} ;

        double Regular(const double x, const int nf) const ;

        double Singular(const double x, const int nf) const {return 0.;} ;
        double Local(const int nf) const {return 0.;} ;
        double SingularIntegrated(const double x, const int nf) const {return 0.;};

    private:
        char entry3_;

        //==========================================================================================//
        //  Convolution between the splitting functions Pgg0/Pqq0
        //  and Pgq0
        //------------------------------------------------------------------------------------------//

        double Pgg0_x_Pgq0(const double x, const int nf) const ;
        double Pqq0_x_Pgq0(const double x) const ;

        //==========================================================================================//
        //  Convolution between the splitting functions Pgq0 and Pqg0
        //------------------------------------------------------------------------------------------//

        double Pgq0_x_Pqg0(const double x, const int nf) const ;
};

class Delta : Function {
    public:
        double Regular(const double x, const int nf) const {return 0.;};
        double Singular(const double x, const int nf) const {return 0.;};
        double Local(const int nf) const {return 1.;};
        double SingularIntegrated(const double x, const int nf) const {return 0.;};
}

#endif
