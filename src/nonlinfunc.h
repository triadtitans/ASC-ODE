#ifndef NONLINFUNC_H
#define NONLINFUNC_H

#include <vector.hpp>
#include <matrix.hpp>


namespace ASC_ode
{
  using namespace ngbla;

  class NonlinearFunction
  {
  public:
    virtual ~NonlinearFunction() = default;
    virtual size_t DimX() const = 0;
    virtual size_t DimF() const = 0;
    virtual void Evaluate (VectorView<double> x, VectorView<double> f) const = 0;
    virtual void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const = 0;
  };


  class IdenticFunction : public NonlinearFunction
  {
    size_t n_;
  public:
    IdenticFunction (size_t n) : n_(n) { } 
    size_t DimX() const override { return n_; }
    size_t DimF() const override { return n_; }
    void Evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      f = x;
    }
    void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
    {
      df = 0.0;
      df.Diag() = 1.0;
    }
  };



  class ConstantFunction : public NonlinearFunction
  {
    Vector<> val_;
  public:
    ConstantFunction (VectorView<double> val) : val_(val) { }
    void Set(VectorView<double> val) { val_ = val; }
    VectorView<double> Get() const { return val_.View(); }
    size_t DimX() const override { return val_.Size(); }
    size_t DimF() const override { return val_.Size(); }
    void Evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      f = val_;
    }
    void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
    {
      df = 0.0;
    }
  };

  
  
  class SumFunction : public NonlinearFunction
  {
    shared_ptr<NonlinearFunction> fa_, fb_;
    double faca_, facb_;
  public:
    SumFunction (shared_ptr<NonlinearFunction> fa,
                 shared_ptr<NonlinearFunction> fb,
                 double faca, double facb)
      : fa_(fa), fb_(fb), faca_(faca), facb_(facb) { } 
    
    size_t DimX() const override { return fa_->DimX(); }
    size_t DimF() const override { return fa_->DimF(); }
    void Evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      fa_->Evaluate(x, f);
      f *= faca_;
      Vector<> tmp(DimF());
      fb_->Evaluate(x, tmp);
      f += facb_*tmp;
    }
    void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
    {
      fa_->EvaluateDeriv(x, df);
      Matrix<> tmp(DimF(), DimX());
      tmp *= faca_;
      fb_->EvaluateDeriv(x, tmp);
      df += facb_*tmp;
    }
  };


  inline auto operator- (shared_ptr<NonlinearFunction> fa, shared_ptr<NonlinearFunction> fb)
  {
    return make_shared<SumFunction>(fa, fb, 1, -1);
  }

  inline auto operator+ (shared_ptr<NonlinearFunction> fa, shared_ptr<NonlinearFunction> fb)
  {
    return make_shared<SumFunction>(fa, fb, 1, 1);
  }

  
  class ScaleFunction : public NonlinearFunction
  {
    shared_ptr<NonlinearFunction> fa_;
    double fac_;
  public:
    ScaleFunction (shared_ptr<NonlinearFunction> fa,
                    double fac)
      : fa_(fa), fac_(fac) { } 
    
    size_t DimX() const override { return fa_->DimX(); }
    size_t DimF() const override { return fa_->DimF(); }
    void Evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      fa_->Evaluate(x, f);
      f *= fac_;

    }
    void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
    {
      fa_->EvaluateDeriv(x, df);
      df *= fac_;
    }
  };


  inline auto operator* (double a, shared_ptr<NonlinearFunction> f)
  {
    return make_shared<ScaleFunction>(f, a);
  }




  // fa(fb)
  class ComposeFunction : public NonlinearFunction
  {
    shared_ptr<NonlinearFunction> fa_, fb_;
  public:
    ComposeFunction (shared_ptr<NonlinearFunction> fa,
                     shared_ptr<NonlinearFunction> fb)
      : fa_(fa), fb_(fb) { } 
    
    size_t DimX() const override { return fb_->DimX(); }
    size_t DimF() const override { return fa_->DimF(); }
    void Evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      Vector<> tmp(fb_->DimF());
      fb_->Evaluate (x, tmp);
      fa_->Evaluate (tmp, f);
    }
    void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
    {
      Vector<> tmp(fb_->DimF());
      fb_->Evaluate (x, tmp);
      
      Matrix<> jaca(fa_->DimF(), fa_->DimX());
      Matrix<> jacb(fb_->DimF(), fb_->DimX());
      
      fb_->EvaluateDeriv(x, jacb);
      fa_->EvaluateDeriv(tmp, jaca);

      df = jaca*jacb;
    }
  };
  
  
  inline auto Compose (shared_ptr<NonlinearFunction> fa, shared_ptr<NonlinearFunction> fb)
  {
    return make_shared<ComposeFunction> (fa, fb);
  }
  
  class EmbedFunction : public NonlinearFunction
  {
    shared_ptr<NonlinearFunction> fa_;
    size_t firstx_, dimx_, firstf_, dimf_;
    size_t nextx_, nextf_;
  public:
    EmbedFunction (shared_ptr<NonlinearFunction> fa,
                   size_t firstx, size_t dimx,
                   size_t firstf, size_t dimf)
      : fa_(fa),
        firstx_(firstx), dimx_(dimx), firstf_(firstf), dimf_(dimf),
        nextx_(firstx+fa->DimX()), nextf_(firstf+fa->DimF())
    { }
    
    size_t DimX() const override { return dimx_; }
    size_t DimF() const override { return dimf_; }
    void Evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      f = 0.0;
      fa_->Evaluate(x.Range(firstx_, nextx_), f.Range(firstf_, nextf_));
    }
    void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
    {
      df = 0;
      fa_->EvaluateDeriv(x.Range(firstx_, nextx_),
                         df.Rows(firstf_, nextf_).Cols(firstx_, nextx_));      
    }
  };

  
  class Projector : public NonlinearFunction
  {
    size_t size_, first_, next_;
  public:
    Projector (size_t size, 
               size_t first, size_t next)
      : size_(size), first_(first), next_(next) { }
    
    size_t DimX() const override { return size_; }
    size_t DimF() const override { return size_; }
    void Evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      f = 0.0;
      f.Range(first_, next_) = x.Range(first_, next_);
    }
    void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
    {
      df = 0.0;
      df.Diag().Range(first_, next_) = 1;
    }
  };

  
}

#endif
