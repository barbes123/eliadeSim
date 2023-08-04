#include "AppsHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include <iomanip>
#include <iostream> 

using namespace std;

G4ThreadLocal G4Allocator<AppsHit>* AppsHitAllocator=0;


AppsHit::AppsHit()
 : G4VHit(),
   fTrackID(-1),
   fChamberNb(-1),
   fEdep(0.),
   fPos(G4ThreeVector())  
{

}


AppsHit::~AppsHit() {}


AppsHit::AppsHit(const AppsHit& right)
  : G4VHit()
{
  fTrackID   = right.fTrackID;
  fChamberNb = right.fChamberNb;
  fEdep      = right.fEdep;
  fPos       = right.fPos;
  fGlobalTime = right.fGlobalTime;
  fFaceEnergy = right.fFaceEnergy;
}


const AppsHit& AppsHit::operator=(const AppsHit& right)
{
  fTrackID   = right.fTrackID;
  fChamberNb = right.fChamberNb;
  fEdep      = right.fEdep;
  fPos       = right.fPos;
  fGlobalTime = right.fGlobalTime;
  fFaceEnergy = right.fFaceEnergy;

  return *this;
}


G4int AppsHit::operator==(const AppsHit& right) const
{
  return ( this == &right ) ? 1 : 0;
}


void AppsHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(fPos);
    circle.SetScreenSize(4.);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}


