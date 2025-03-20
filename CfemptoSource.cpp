#include "CfemptoSource.h"
#include "TRandom.h"

float femptoSource::getKstarFinal(float coalProb) const {
  float Efin = TMath::Max(0., (kineticSource() + potentialSource() - EBOUND*coalProb)/(1 - coalProb));

  return sqrt(Efin*2*MRED);
}
//_________________________________________________________________________
void femptoSource::setCharges(float cS, float cC){
  mSourceV->SetParameter(2,cS);
  mSourceV->SetParameter(3,cC);

  mDeuteronV->SetParameter(5,cS);
  mDeuteronV->SetParameter(6,cC);
}
//_________________________________________________________________________
void femptoSource::setKstar(float kstar) {
  if(! mIsInitialized){
    init();
  }
  if(mSourceRadius < 0){
    static const float factor = sqrt(3*0.5) * HCUT;
    float radius = factor/kstar;
    setSourceRadius(-radius);
    mSourceKin->SetParameter(2, 0);
    mCoalescenceRe->SetParameter(2, 0);
    mCoalescenceIm->SetParameter(2, 0);
  } else {
    mSourceKin->SetParameter(2, kstar);
    mCoalescenceRe->SetParameter(2, kstar);
    mCoalescenceIm->SetParameter(2, kstar);
    setSourceRadius(mSourceRadius); // to trigger new normalization
  }

}
//_________________________________________________________________________
void femptoSource::setSourceRadius(float radius) {
  if(! mIsInitialized){
    init();
  }
  if(radius > 0){
    mSourceRadius = radius;
  } else {
    mSourceRadius = -1;
    radius *= -1;
  }

  if(radius > 0){
    mSource2int->SetParameter(1,radius);
    mSource->SetParameter(1,radius);
    mUSource->SetParameter(1,radius);
    mSourceKin->SetParameter(1,radius);
    mSourceV->SetParameter(1,radius);
    mCoalescenceRe->SetParameter(1,radius);
    mCoalescenceIm->SetParameter(1,radius);
  } else {
    mSource2int->SetParameter(1,0);
    mSource->SetParameter(1,0);
    mUSource->SetParameter(1,0);
    mSourceKin->SetParameter(1,0);
    mSourceV->SetParameter(1,0);
    mCoalescenceRe->SetParameter(1,0);
    mCoalescenceIm->SetParameter(1,0);
  }

  mSource2int->SetParameter(0,1);

  mMaxIntRange = TMath::Max(20., double(radius*3));
  if(mMaxIntRange > 100){
    mMaxIntRange = 100;
  }

  double err;
  double norm = 1./sqrt(mSource2int->IntegralOneDim(0,mMaxIntRange,1E-8,1E-8,err));

  mSource->SetParameter(0,norm);
  mUSource->SetParameter(0,norm);
  mSource2int->SetParameter(0,norm);
  mSourceKin->SetParameter(0,norm);
  mSourceV->SetParameter(0,norm);
  mCoalescenceRe->SetParameter(0,norm);
  mCoalescenceIm->SetParameter(0,norm);
}
//_________________________________________________________________________
double femptoSource::coalescenceRe(double *x,double *pm){
  double r = x[0];
  double pDeu[5] = {pm[3], pm[4], pm[5], pm[6], pm[7]};

  return waveDeuteron(x,pDeu) * sourceCos(x,pm) * 4 * TMath::Pi() * r * r;
}
//_________________________________________________________________________
double femptoSource::coalescenceIm(double *x,double *pm){
  double r = x[0];
  double pDeu[5] = {pm[3], pm[4], pm[5], pm[6], pm[7]};

  return waveDeuteron(x,pm) * sourceSin(x,pm) * 4 * TMath::Pi() * r * r;
}
//_________________________________________________________________________
double femptoSource::uDeuteron(double *x,double *pm){
    double r = x[0];

    double res = pm[0];
    double radius = pm[1];
    double k1 = pm[2];
    double k2 = pm[3];
    double normRight = pm[4];

    if(r < radius){
      res *= sin(k1*r);
    } else {
      res *= TMath::Exp(-k2*r) * normRight;
    }

    return res;
}
//_________________________________________________________________________
double femptoSource::waveDeuteron(double *x,double *pm){
    double r = x[0];

    return uDeuteron(x,pm) / r;
}
//_________________________________________________________________________
double femptoSource::intDeuteron(double *x,double *pm){
    double r = x[0];

    return waveDeuteron(x,pm) * waveDeuteron(x,pm) * 4 * TMath::Pi() * r * r;
}
//_________________________________________________________________________
double femptoSource::deuteronKin(double *x,double *pm){
    double r = x[0];

    double res = 0;
    double radius = pm[1];
    double k1 = pm[2];
    double k2 = pm[3];

    if(r < radius){
      res = k1*k1*HCUT*HCUT*0.5/MRED;
    } else {
      res = -k2*k2*HCUT*HCUT*0.5/MRED;
    }

  return res * intDeuteron(x,pm);
}
//_________________________________________________________________________
double femptoSource::deuteronV(double *x,double *pm){
  double r = x[0];
  float chargeStrong = pm[5];
  float chargeColoumb = pm[6];
  double val = 0;
  if(r < mStrongR){
    val = - mStrong * chargeStrong;
  }
  val += chargeColoumb * mCoulomb / r;

  return val * intDeuteron(x,pm);
}
//_________________________________________________________________________
double femptoSource::source(double *x,double *pm){
    double r = x[0];
    double norm = pm[0];
    double radius = pm[1];
    return norm * TMath::Exp(-r*r*0.5/(radius*radius));
}
//_________________________________________________________________________
double femptoSource::sourceCos(double *x,double *pm){
    double r = x[0];
    double norm = pm[0];
    double radius = pm[1];
    double kstar = pm[2];

    return norm * TMath::Exp(-r*r*0.5/(radius*radius)) * cos(kstar*r/HCUT);
}
//_________________________________________________________________________
double femptoSource::sourceSin(double *x,double *pm){
    double r = x[0];
    double norm = pm[0];
    double radius = pm[1];
    double kstar = pm[2];

    return norm * TMath::Exp(-r*r*0.5/(radius*radius)) * sin(kstar*r/HCUT);
}
//_________________________________________________________________________
double femptoSource::usource(double *x,double *pm){
    double r = x[0];

    return source(x,pm) * r;
}
//_________________________________________________________________________
double femptoSource::source2int(double *x,double *pm){
    double r = x[0];

    return source(x,pm) * source(x,pm) * 4 * TMath::Pi() * r * r;
}
//_________________________________________________________________________
double femptoSource::sourceKin(double *x,double *pm){
  double r = x[0];
  double radius = pm[1];
  double kstar = pm[2];
  double R02inv = 1./(radius*radius);

  double val = 3 * R02inv * (HCUT*HCUT)/(2*MRED);
  val -= r*r*R02inv*R02inv * (HCUT*HCUT)/(2*MRED);
  val += kstar*kstar / (2*MRED);

  return val * source2int(x,pm);
}
//_________________________________________________________________________
double femptoSource::sourceV(double *x,double *pm){
  double r = x[0];
  float chargeStrong = pm[2];
  float chargeColoumb = pm[3];
  double val = - (r < mStrongR) * mStrong * chargeStrong;
  val += chargeColoumb * mCoulomb / r;

  return val * source2int(x,pm);
}
//_________________________________________________________________________
void femptoSource::init(){
  mK1 = sqrt(2*MRED*(EBOUND+mStrong))/HCUT;
  mK2 = sqrt(-2*MRED*EBOUND)/HCUT;
  printf("Deuteron k1=%f - k2=%f - Radius=%f - V = %f (E=%f)\n",mK1,mK2,mStrongR,mStrong,EBOUND);
  mNormRight = sin(mK1*mStrongR) * TMath::Exp(mK2*mStrongR);

  mSource = new TF1("fSource",femptoSource::source,0,20,2);
  mUSource = new TF1("fUSource",femptoSource::usource,0,20,2);
  mSource2int = new TF1("fSource2int",femptoSource::source2int,0,20,2);
  mSourceKin = new TF1("fSourceKin",femptoSource::sourceKin,0,20,3);
  mSourceV = new TF1("fSourceV",femptoSource::sourceV,0,20,4);

  mDeuteron2int = new TF1("fDeuteron2int",femptoSource::intDeuteron,0,20,5);
  mDeuteron2int->SetParameter(0,1);
  mDeuteron2int->SetParameter(1,mStrongR);
  mDeuteron2int->SetParameter(2,mK1);
  mDeuteron2int->SetParameter(3,mK2);
  mDeuteron2int->SetParameter(4,mNormRight);
  double err;
  double norm = 1./sqrt(mDeuteron2int->IntegralOneDim(0,20,1E-8,1E-8,err));
  mDeuteron2int->SetParameter(0,norm);

  mDeuteron = new TF1("fDeuteron",femptoSource::waveDeuteron,0,20,5);
  mDeuteron->SetParameter(0,norm);
  mDeuteron->SetParameter(1,mStrongR);
  mDeuteron->SetParameter(2,mK1);
  mDeuteron->SetParameter(3,mK2);
  mDeuteron->SetParameter(4,mNormRight);

  mUDeuteron = new TF1("fUDeuteron",femptoSource::uDeuteron,0,20,5);
  mUDeuteron->SetParameter(0,norm);
  mUDeuteron->SetParameter(1,mStrongR);
  mUDeuteron->SetParameter(2,mK1);
  mUDeuteron->SetParameter(3,mK2);
  mUDeuteron->SetParameter(4,mNormRight);

  mDeuteronKin = new TF1("fDeuteronKin",femptoSource::deuteronKin,0,20,5);
  mDeuteronKin->SetParameter(0,norm);
  mDeuteronKin->SetParameter(1,mStrongR);
  mDeuteronKin->SetParameter(2,mK1);
  mDeuteronKin->SetParameter(3,mK2);
  mDeuteronKin->SetParameter(4,mNormRight);

  mDeuteronV = new TF1("fDeuteronV",femptoSource::deuteronV,0,20,7);
  mDeuteronV->SetParameter(0,norm);
  mDeuteronV->SetParameter(1,mStrongR);
  mDeuteronV->SetParameter(2,mK1);
  mDeuteronV->SetParameter(3,mK2);
  mDeuteronV->SetParameter(4,mNormRight);

  mCoalescenceRe = new TF1("fCoalescenceRe",femptoSource::coalescenceRe,0,20,8);
  mCoalescenceRe->SetParameter(3, norm);
  mCoalescenceRe->SetParameter(4, mStrongR);
  mCoalescenceRe->SetParameter(5, mK1);
  mCoalescenceRe->SetParameter(6, mK2);
  mCoalescenceRe->SetParameter(7, mNormRight);
  mCoalescenceIm = new TF1("fCoalescenceIm",femptoSource::coalescenceIm,0,20,8);
  mCoalescenceIm->SetParameter(3, norm);
  mCoalescenceIm->SetParameter(4, mStrongR);
  mCoalescenceIm->SetParameter(5, mK1);
  mCoalescenceIm->SetParameter(6, mK2);
  mCoalescenceIm->SetParameter(7, mNormRight);

  vfempto::init();

  if(mSourceRadius >= 0){
    setSourceRadius(mSourceRadius);
  } else {
    setSourceRadius(0);
  }
}
//_________________________________________________________________________
float femptoSource::calcProb(){
  double err;
  float probRe = mCoalescenceRe->IntegralOneDim(0,mMaxIntRange,1E-8,1E-8,err);
  probRe *= probRe;
  float probIm = mCoalescenceIm->IntegralOneDim(0,mMaxIntRange,1E-8,1E-8,err);
  probIm *= probIm;

  return (probRe + probIm) * mSpinCoalFactor;
}
//_________________________________________________________________________
float femptoSource::getCoalProb(const particleMC& p1, const particleMC& p2) {
  double kstar = utils::getKstar(p1,p2) * 1E3;
  setKstar(kstar);

  return calcProb();
}
//_________________________________________________________________________
double femptoSource::doInteract(particleMC& p1, particleMC& p2, float chargeColoumb, float chargeStrong, float sumRadii, float *pos, float *posLab){
  double kstar = utils::getKstar(p1,p2);
  setKstar(kstar);

  TLorentzVector pSum = p1.q + p2.q;
  TVector3 b = pSum.BoostVector();
  TVector3 bInv = -b;

  p1.q.Boost(bInv);
  p2.q.Boost(bInv);

  setCharges(chargeStrong,chargeColoumb);

  float coalProb = 0;
  if(std::abs(chargeStrong) > 0.9 && std::abs(chargeColoumb) < 1E-3) { // coalescence allowed
    coalProb = calcProb();
  }

  double momFinal = getKstarFinal(coalProb);
  float scaling = momFinal / p1.q.P();

  double m1 = p1.q.M();
  double px1= p1.q.Px()*scaling;
  double py1= p1.q.Py()*scaling;
  double pz1= p1.q.Pz()*scaling;
  double e1 = sqrt(px1*px1 + py1*py1 + pz1*pz1 + m1*m1);
  double ptEx = std::abs(momFinal) - std::abs(p1.q.P());
  p1.q.SetPxPyPzE(px1,py1,pz1,e1);
  double m2 = p2.q.M();
  double px2= p2.q.Px()*scaling;
  double py2= p2.q.Py()*scaling;
  double pz2= p2.q.Pz()*scaling;
  double e2 = sqrt(px2*px2 + py2*py2 + pz2*pz2 + m2*m2);
  p2.q.SetPxPyPzE(px2,py2,pz2,e2);

  p1.q.Boost(b);
  p2.q.Boost(b);

  return ptEx;
}
