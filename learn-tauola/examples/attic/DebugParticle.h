#ifndef PART_TESTING_ROUTINES
#define PART_TESTING_ROUTINES
#include <iostream>
#include <math.h>
using std::ostream;
using std::cout;
using std::endl;
using std::ios;

class DebugParticle
{
public:
	DebugParticle():Px(0.0),Py(0.0),Pz(0.0),E(0.0),PdgId(0) {}
	DebugParticle(const DebugParticle& x):Px(x.Px),Py(x.Py),Pz(x.Pz),E(x.E),PdgId(x.PdgId) {}
	DebugParticle(double x,double y, double z, double e):Px(x),Py(y),Pz(z),E(e),PdgId(0) {}
	DebugParticle(double x,double y, double z, double e,int id):Px(x),Py(y),Pz(z),E(e),PdgId(id) {}
public:
	const DebugParticle& set(double x, double y, double z, double e)
	{
		Px=x;
		Py=y;
		Pz=z;
		E=e;
		return *this;
	}
	const bool isZero() const
	{
		return ( (*this)*(*this) < 0.00000001 ) ? true : false ;
	}
	const double GetM() const
	{
		double sq = E*E-Px*Px-Py*Py-Pz*Pz;
		if(sq<0) return -sqrt(-sq);
		return sqrt(sq);
	}
public:
	const bool operator==(const DebugParticle& x) const
	{
		DebugParticle buf=*this;
		buf-=x;
		return (buf*buf<0.00000001) ? true : false ;
	}
	const bool operator!=(const DebugParticle& x) const
	{
		return !(*this==x);
	}
	const DebugParticle& operator=(const DebugParticle& x)
	{
		Px=x.Px;
		Py=x.Py;
		Pz=x.Pz;
		E=x.E;
		PdgId=x.PdgId;
		return *this;
	}
	const DebugParticle& operator+=(const DebugParticle& x)
	{
		Px+=x.Px;
		Py+=x.Py;
		Pz+=x.Pz;
		E+=x.E;
		return *this;
	}
	const DebugParticle& operator-=(const DebugParticle& x)
	{
		Px-=x.Px;
		Py-=x.Py;
		Pz-=x.Pz;
		E-=x.E;
		return *this;
	}
	const DebugParticle operator+(const DebugParticle &x) const
	{
		DebugParticle buf(*this);
		buf.PdgId=0;
		return buf+=x;
	}
	const DebugParticle operator-(const DebugParticle &x) const
	{
		DebugParticle buf(*this);
		buf.PdgId=0;
		return buf-=x;
	}
	const double operator*(const DebugParticle &x)   const
	{
		return dot_product(x);
	}
	const double operator|(const DebugParticle &x)   const
	{
		return virtuality(x);
	}
	const DebugParticle& operator>>(const DebugParticle &x)
	{
		return Boost(x);
	}
public:
	const double dot_product(const DebugParticle &x) const
	{
		return Px*x.Px+Py*x.Py+Pz*x.Pz;
	}
	const double virtuality(const DebugParticle &x)  const
	{
		DebugParticle diff=*this-x;
		return diff*diff-(*this)*(*this)-x*x;
	}
	const double length() const
	{
		return sqrt(Px*Px+Py*Py+Pz*Pz);
	}
	const DebugParticle& Boost(const DebugParticle& x)
	{
		return Boost(x.Px,x.Py,x.Pz,x.E,x.GetM());
	}
	const DebugParticle& Boost(double x, double y, double z,double e, double m)
	{
		double betx=-x/m;
		double bety=-y/m;
		double betz=-z/m;
		double gam=e/m;
		double pb=betx*Px+bety*Py+betz*Pz;
		Px=Px+betx*(E+pb/(gam+1.0));
		Py=Py+bety*(E+pb/(gam+1.0));
		Pz=Pz+betz*(E+pb/(gam+1.0));
		E=E*gam+pb;
		return *this;
	}
	const DebugParticle& smallestVirtuality(const DebugParticle &p1, const DebugParticle &p2, const DebugParticle &q1, const DebugParticle &q2) const
	{
		const DebugParticle *ret=&p1;
		double virt=*this|p1;
		double subVirt=*this|p2;
		if( subVirt < virt ) { ret=&p2; virt=subVirt; }
		subVirt=*this|q1;
		if( subVirt < virt ) { ret=&q1; virt=subVirt; }
		subVirt=*this|q2;
		if( subVirt < virt ) { ret=&q2; virt=subVirt; }
		return *ret;
	}
public:
	static const double S(DebugParticle &p1, DebugParticle &p2, DebugParticle &q1, DebugParticle &q2)
	{
		DebugParticle pp=p1+p2;
		DebugParticle qq=q1+q1;
		if(pp*pp!=qq*qq) cout<<"DebugParticle::S - Not equal: "<<pp*pp<<" "<<qq*qq<<endl;
		return pp*pp;
	}
	const double cosTheta(DebugParticle &t)
	{
		double cT = (*this)*t/(length()*t.length());
		if(t.PdgId>0) return -cT;
		return cT;
	}
public:
	double Px,Py,Pz,E;
	int PdgId;
};

ostream& operator<<(ostream &s, const DebugParticle& p)
{
	int lastP=s.precision(3);
	s.setf(ios::showpos|ios::scientific);
	s<<p.Px<<",";
	s.width(14);
	s<<p.Py<<",";
	s.width(14);
	s<<p.Pz<<",";
	s.width(14);
	s<<p.E<<",";
	s.width(14);
	s<<p.GetM();
	if(p.PdgId) s<<"  ("<<p.PdgId<<")"<<endl;
	else s<<endl;
	s.precision(lastP);
	s.unsetf(ios::showpos|ios::scientific);
	return s;
}

#endif
/*
int main()  //self-test
{
	DebugParticle x;
	DebugParticle y(1,2,3,4);
	DebugParticle z(2,3,4,8,-112);
	DebugParticle w=z;
	cout<<x<<y<<z<<w<<endl;
	cout<<y-z<<y<<endl;
	z-=y;
	y+=z;
	cout<<y.smallestVirtuality(x,y+x,z,w)<<endl;
	cout<<"Y*Y: "<<y*y<<endl;                     //dot_product
	cout<<"X|Y: "<<(x|y)<<endl;                   //virtuality
	cout<<(x>>y)<<endl;                           //boost
}
*/

/*
//test inside select()
    DebugParticle x(tauPlus->GetPx(),tauPlus->GetPy(),tauPlus->GetPz(),tauPlus->GetE(),tauPlus->GetPDGId());
    DebugParticle y(tauMinus->GetPx(),tauMinus->GetPy(),tauMinus->GetPz(),tauMinus->GetE(),tauMinus->GetPDGId());
    DebugParticle z(grand1->GetPx(),grand1->GetPy(),grand1->GetPz(),grand1->GetE(),grand1->GetPDGId());
    DebugParticle w(grand2->GetPx(),grand2->GetPy(),grand2->GetPz(),grand2->GetE(),grand2->GetPDGId());
    cout<<x<<y<<z<<w<<endl;
    cout<<"Czy sie zgadza? "<<(x+y==w+z)<<endl;
    //boost x to x+y;
    DebugParticle buf=x+y;
    //Like that:
    x.Boost(buf);
    y.Boost(buf);
    //or
    z>>buf;
    w>>buf;
*/
