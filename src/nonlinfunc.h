#ifndef NONLINFUNC_H
#define NONLINFUNC_H

#include "vector.h"
#include "matrix.h"


namespace ASC_ode
{
  using namespace ASC_bla;

  class NonlinearFunction
  {
  public:
    virtual ~NonlinearFunction() = default;
    virtual size_t DimX() const = 0;
    virtual size_t DimF() const = 0;
    virtual void Evaluate (VectorView<double> x, VectorView<double> f) const = 0;
    virtual void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const = 0;
  };

  class IdentityFunction : public NonlinearFunction
  {
    size_t n;
  public:
    IdentityFunction (size_t _n) : n(_n) { } 
    size_t DimX() const override { return n; }
    size_t DimF() const override { return n; }
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

  class LinearFunction : public NonlinearFunction
  {

    Matrix<double> m_;
    public:
    template <typename T>
    LinearFunction (const MatrixExpr<T>& m) : m_(m) { 
    } 
    size_t DimX() const override { return m_.Width(); }
    size_t DimF() const override { return m_.Height(); }
    void Evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      f = m_*x;
    }
    
    void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
    {
      df = m_;
    }
  };

  class StackedFunction : public NonlinearFunction
  {
    std::vector<std::shared_ptr<NonlinearFunction>> _functions;
    public:
    StackedFunction(){};
    void addFunction(std::shared_ptr<NonlinearFunction> func){_functions.push_back(func);}
    size_t DimX() const override {
      size_t sum=0;
      for(auto ptr : _functions)
        sum += ptr->DimX();
      return sum;
    }
    size_t DimF() const override {
      size_t sum=0;
      for(auto ptr : _functions)
        sum += ptr->DimF();
      return sum;
    }
    void Evaluate (VectorView<double> x, VectorView<double> f) const override{
      size_t cursor_x=0;
      size_t cursor_f=0;
      for(auto func : _functions){
        func->Evaluate(x.Range(cursor_x,func->DimX()),f.Range(cursor_f,func->DimF()));
        cursor_f+=func->DimF();
        cursor_x+=func->DimX();
      }
    }
    void EvaluateDeriv (VectorView<double> x, MatrixView<double> f) const override{
      size_t cursor_x=0;
      size_t cursor_f=0;
      f=0;
      for(auto func : _functions){
        //The Jacobian of the stacked function is a block diagonal matrix, consisting of the individual Jacobians
        MatrixView<double> currentBlock = f.Cols(cursor_x,func->DimX()).Rows(cursor_f,func->DimF());
        func->EvaluateDeriv(x.Range(cursor_x,func->DimX()),currentBlock);
        cursor_f+=func->DimF();
        cursor_x+=func->DimX();
      }
    }
  };

  void dNumeric(const NonlinearFunction& f, VectorView<double> x, MatrixView<double> df){
    double eps = 1e-8;
    Vector<> xl(f.DimX()), xr(f.DimX()), fl(f.DimF()), fr(f.DimF());
    for (size_t i = 0; i < f.DimX(); i++)
      {
        xl = x;
        xl(i) -= eps;
        xr = x;
        xr(i) += eps;
        f.Evaluate (xl, fl);
        f.Evaluate (xr, fr);
        df.Col(i) = 1/(2*eps) * (fr-fl);
      }
  }

class NumericDerivative : public NonlinearFunction
{
  std::shared_ptr<NonlinearFunction> g_;
  public:
  size_t DimX() const override { return 18; }
  size_t DimF() const override { return 18; }
  NumericDerivative(std::shared_ptr<NonlinearFunction> g): g_(g){};
  void Evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    Matrix<double> fmat (1, 18);
    g_->EvaluateDeriv(x, fmat);
    f = fmat.Row(0);
  }
  virtual void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const
  {
    dNumeric(*this,x,df);
  }
};

  class ConstantFunction : public NonlinearFunction
  {
    Vector<> val;
  public:
    ConstantFunction (VectorView<double> _val) : val(_val) { }
    void Set(VectorView<double> _val) { val = _val; }
    VectorView<double> Get() const { return val.View(); }
    size_t DimX() const override { return val.Size(); }
    size_t DimF() const override { return val.Size(); }
    void Evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      f = val;
    }
    void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
    {
      df = 0.0;
    }
  };

  class SumFunction : public NonlinearFunction
  {
    std::shared_ptr<NonlinearFunction> fa, fb;
    double faca, facb;
  public:
    SumFunction (std::shared_ptr<NonlinearFunction> _fa,
                 std::shared_ptr<NonlinearFunction> _fb,
                 double _faca, double _facb)
      : fa(_fa), fb(_fb), faca(_faca), facb(_facb) { } 
    
    size_t DimX() const override { return fa->DimX(); }
    size_t DimF() const override { return fa->DimF(); }
    void Evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      fa->Evaluate(x, f);
      f *= faca;
      Vector<> tmp(DimF());
      fb->Evaluate(x, tmp);
      f += facb*tmp;
    }
    void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
    {
      fa->EvaluateDeriv(x, df);
      df *= faca;
      Matrix<> tmp(DimF(), DimX());
      fb->EvaluateDeriv(x, tmp);
      df += facb*tmp;
    }
  };

  inline auto operator- (std::shared_ptr<NonlinearFunction> fa, std::shared_ptr<NonlinearFunction> fb)
  {
    return std::make_shared<SumFunction>(fa, fb, 1, -1);
  }

  inline auto operator+ (std::shared_ptr<NonlinearFunction> fa, std::shared_ptr<NonlinearFunction> fb)
  {
    return std::make_shared<SumFunction>(fa, fb, 1, 1);
  }

  class ScaleFunction : public NonlinearFunction
  {
    std::shared_ptr<NonlinearFunction> fa;
    double fac;
  public:
    ScaleFunction (std::shared_ptr<NonlinearFunction> _fa,
                   double _fac)
      : fa(_fa), fac(_fac) { } 
    
    size_t DimX() const override { return fa->DimX(); }
    size_t DimF() const override { return fa->DimF(); }
    void Evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      fa->Evaluate(x, f);
      f *= fac;

    }
    void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
    {
      fa->EvaluateDeriv(x, df);
      df *= fac;
    }
  };

  inline auto operator* (double a, std::shared_ptr<NonlinearFunction> f)
  {
    return std::make_shared<ScaleFunction>(f, a);
  }

  // fa(fb)
  class ComposeFunction : public NonlinearFunction
  {
    std::shared_ptr<NonlinearFunction> fa, fb;
  public:
    ComposeFunction (std::shared_ptr<NonlinearFunction> _fa,
                     std::shared_ptr<NonlinearFunction> _fb)
      : fa(_fa), fb(_fb) { } 
    
    size_t DimX() const override { return fb->DimX(); }
    size_t DimF() const override { return fa->DimF(); }
    void Evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      Vector<> tmp(fb->DimF());
      fb->Evaluate (x, tmp);
      fa->Evaluate (tmp, f);
    }
    void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
    {
      Vector<> tmp(fb->DimF());
      fb->Evaluate (x, tmp);
      
      Matrix<> jaca(fa->DimF(), fa->DimX());
      Matrix<> jacb(fb->DimF(), fb->DimX());
      
      fb->EvaluateDeriv(x, jacb);
      fa->EvaluateDeriv(tmp, jaca);

      df = jaca*jacb;
    }
  };
  
  
  inline auto Compose (std::shared_ptr<NonlinearFunction> fa, std::shared_ptr<NonlinearFunction> fb)
  {
    return std::make_shared<ComposeFunction> (fa, fb);
  }
  
  class EmbedFunction : public NonlinearFunction
  {
    std::shared_ptr<NonlinearFunction> fa;
    size_t firstx, dimx, firstf, dimf;
    size_t nextx, nextf;
  public:
    EmbedFunction (std::shared_ptr<NonlinearFunction> _fa,
                   size_t _firstx, size_t _dimx,
                   size_t _firstf, size_t _dimf)
      : fa(_fa),
        firstx(_firstx), dimx(_dimx), firstf(_firstf), dimf(_dimf),
        nextx(_firstx+_fa->DimX()), nextf(_firstf+_fa->DimF())
    { }
    
    size_t DimX() const override { return dimx; }
    size_t DimF() const override { return dimf; }
    void Evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      f = 0.0;
      fa->Evaluate(x.Range(firstx, nextx), f.Range(firstf, nextf));
    }
    void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
    {
      df = 0;
      fa->EvaluateDeriv(x.Range(firstx, nextx),
                        df.Rows(firstf, nextf).Cols(firstx, nextx));      
    }
  };

  
  class Projector : public NonlinearFunction
  {
    size_t size, first, next;
  public:
    Projector (size_t _size, 
               size_t _first, size_t _next)
      : size(_size), first(_first), next(_next) { }
    
    size_t DimX() const override { return size; }
    size_t DimF() const override { return size; }
    void Evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      f = 0.0;
      f.Range(first, next) = x.Range(first, next);
    }
    void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
    {
      df = 0.0;
      df.Diag().Range(first, next) = 1;
    }
  };

  class BlockFunction : public NonlinearFunction
  {
    std::shared_ptr<NonlinearFunction> func;
    size_t s;
  public:
    BlockFunction (size_t s, size_t n) : s(s) { }

    size_t DimX() const override { return s*func->DimX(); }
    size_t DimF() const override { return s*func->DimF(); }

    void Evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      size_t fdim = func->DimF();
      size_t xdim = func->DimX();
      for (size_t i = 0; i<s; i++) {
        func->Evaluate(x.Range(i*xdim, (i+1)*xdim), f.Range(i*fdim, (i+1)*fdim));
      }
    }
    void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
    {
      df = 0.0;
      size_t fdim = func->DimF();
      size_t xdim = func->DimX();
      for (size_t i = 0; i<s; i++) {
        func->EvaluateDeriv(x.Range(i*xdim,(i+1)*xdim), df.Rows(i*fdim,fdim).Cols(i*xdim, xdim));
      }
    }
  };
}

#endif
