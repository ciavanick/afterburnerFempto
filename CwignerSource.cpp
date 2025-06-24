#include "CwignerSource.h"
#include "CwignerUtils.h"

void wignerSource::initFunctions()
{
    wignerUtils::init();
}
//_________________________________________________________________________
void wignerSource::setFunctionsParameters()
{
    wignerUtils::setFunctionsParameters();
}
//_________________________________________________________________________
/*
void wignerSource::setRadius(double radius)
{
    mRadius = radius;
    mWxJ->SetParameter(0, 1.);
    reSetRadius();
    normalization();
    reSetNorm();
}
    */
//_________________________________________________________________________
void wignerSource::setR0(double r0)
{
    wignerUtils::setSourceRadius(r0);
}
//_________________________________________________________________________
void wignerSource::setRadiusK(double k)
{
    wignerUtils::setRadiusK(k);
}
//_________________________________________________________________________
// to check if it is useful
void wignerSource::setKstar(double k)
{
    wignerUtils::setKstar(k, 0.);
}
//_________________________________________________________________________
/*
void wignerSource::setKIn(double k)
{
    mKin = k;
}
*/
//_________________________________________________________________________
void wignerSource::setRanges(double xmin, double ymin, double xmax, double ymax)
{
    wignerUtils::setRanges(xmin, ymin, xmax, ymax);
}
//_________________________________________________________________________
void wignerSource::setMu(double mu)
{
    wignerUtils::setMu(mu);
}
//_________________________________________________________________________
void wignerSource::setRWidth(double rWidth)
{
    wignerUtils::setRWidth(rWidth);
}
//_________________________________________________________________________
void wignerSource::setV0(double v0)
{
    wignerUtils::setV0(v0);
}
//_________________________________________________________________________
TF2 *wignerSource::getWignerFunction()
{
    return wignerUtils::getWignerFunction();
}
//_________________________________________________________________________
TF2 *wignerSource::getWignerFunctionForItself()
{
    return wignerUtils::getWignerFunctionForItself();
}
//_________________________________________________________________________
TF2 *wignerSource::getWignerFunctionForJacobian()
{
    return wignerUtils::getWignerFunctionForJacobian();
}
//_________________________________________________________________________
TF2 *wignerSource::getWignerFunction2ForJacobian()
{
    return wignerUtils::getWignerFunction2ForJacobian();
}
//_________________________________________________________________________
TF2 *wignerSource::getKineticEnergyFunction()
{
    return wignerUtils::getKineticEnergyFunction();
}
//_________________________________________________________________________
TF2 *wignerSource::getPotentialEnergyFunction()
{
    return wignerUtils::getPotentialEnergyFunction();
}
//_________________________________________________________________________
TF2 *wignerSource::getHamiltonianFunction()
{
    return wignerUtils::getHamiltonianFunction();
}
//_________________________________________________________________________
TF2 *wignerSource::getWignerKinetic()
{
    return wignerUtils::getWignerKinetic();
}
//_________________________________________________________________________
TF2 *wignerSource::getWignerPotential()
{
    return wignerUtils::getWignerPotential();
}
//_________________________________________________________________________
TF2 *wignerSource::getWignerHamiltonan()
{
    return wignerUtils::getWignerHamiltonan();
}
//_________________________________________________________________________
TF2 *wignerSource::getWignerDeuteron()
{
    return wignerUtils::getWignerDeuteron();
}
//_________________________________________________________________________
TF2 *wignerSource::getWignerDeuteronIntegral()
{
    return wignerUtils::getWignerDeuteronIntegral();
}
//_________________________________________________________________________
TF2 *wignerSource::getCoalescenceProbability()
{
    return wignerUtils::getCoalescenceProbability();
}
//_________________________________________________________________________
double wignerSource::getNorm()
{
    return wignerUtils::getNorm();
}
//_________________________________________________________________________
double wignerSource::getwK()
{
    return wignerUtils::kineticSource();
}
//_________________________________________________________________________
double wignerSource::getwV()
{
    return wignerUtils::potentialSource();
}
//_________________________________________________________________________
double wignerSource::getwH()
{
    return wignerUtils::getWH();
}
//_________________________________________________________________________
double wignerSource::checkWxW()
{
    return wignerUtils::checkWxW();
}
//_________________________________________________________________________
double wignerSource::getRadius()
{
    return wignerUtils::getRadius();
}
//_________________________________________________________________________
double wignerSource::getKStar()
{
    return wignerUtils::getKStar();
}
//_________________________________________________________________________
double wignerSource::getRMin()
{
    return wignerUtils::getRMin();
}
//_________________________________________________________________________
double wignerSource::getRMax()
{
    return wignerUtils::getRMax();
}
//_________________________________________________________________________
double wignerSource::getPMin()
{
    return wignerUtils::getPMin();
}
//_________________________________________________________________________
double wignerSource::getPMax()
{
    return wignerUtils::getPMax();
}
//_________________________________________________________________________
double wignerSource::getMu()
{
    return wignerUtils::getMu();
}
//_________________________________________________________________________
double wignerSource::getRWidth()
{
    return wignerUtils::getRWidth();
}
//_________________________________________________________________________
double wignerSource::getV0()
{
    return wignerUtils::getV0();
}
//_________________________________________________________________________
double wignerSource::getcoal()
{
    return wignerUtils::getcoal();
}
//_________________________________________________________________________
double wignerSource::getDeuteronInt()
{
    return wignerUtils::getDeuteronInt();
}
//_________________________________________________________________________
double wignerSource::doInteract(particleMC& p1, particleMC& p2){
    return 0.;
}
//_________________________________________________________________________
float wignerSource::getCoalProb(const particleMC& p1, const particleMC& p2){
    return 0.;
}
