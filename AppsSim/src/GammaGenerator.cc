#include "GammaGenerator.hh"

// Container class for gammas as needed in PrimaryGeneratorAction
MCGamma::MCGamma() : energy(0. * MeV), time(0. * ns), position(G4ThreeVector(0. * mm, 0. * mm, 0. * mm)), polarization(G4ThreeVector(0., 0., 0.)), direction(G4ParticleMomentum(0. * MeV, 0. * MeV, 0. * MeV)) {}

MCGamma::~MCGamma() {}

GammaGenerator::GammaGenerator(CLHEP::HepRandomEngine *rgen) : randGen(rgen), Electron4Mom(G4LorentzVector()), Photon4Mom(G4LorentzVector()), Gamma4Mom(G4LorentzVector()), Electron3Pos(G4ThreeVector()), Photon3Pol(G4ThreeVector()), GammaStokes(G4ThreeVector()), Electron3Boost(G4ThreeVector())
{
  fGamma = new MCGamma();
  flatRand = new G4RandFlat(randGen);
  expRand = new G4RandExponential(randGen);
  gaussRand = new G4RandGauss(randGen);

  nGamma = 0;
  XSum = 0.;
  mec2 = electron_mass_c2;

  std::ifstream inFile;
  G4String inFileName = "GammaSource.inp";
  inFile.open(inFileName, ios::in | ios::binary);
  if (!inFile)
  {
    G4cout << "GammaGenerator - Error: Can't open the input file." << G4endl;
    exit(1);
  }
  else
  {
    G4cout << "GammaGenerator: reading input from " << inFileName << G4endl;
  }
  inFile >> eEner >> deEner >> eSigZ >> eEmit >> eAlpha >> eBeta >> eNB >> eChar;
  inFile >> lLamb >> dlEner >> lRayl >> lSigT >> lPowr >> lPolD >> lPolA >> elCross >> elFreq;
  inFile.close();
  eEner *= MeV;
  deEner *= perCent;
  eSigZ *= mm;
  eEmit *= nanometer;
  eBeta *= m;
  eChar *= (1.e-09 * coulomb);
  lLamb *= nanometer;
  dlEner *= perCent;
  lRayl *= m;
  lSigT *= picosecond;
  lPowr *= (1.e+03 * watt);
  lPolA *= deg;
  elCross *= mrad;
  elFreq *= megahertz;
  if (lPolD < 0 || lPolD > 1)
  {
    G4cout << "GammaGenerator - ERROR: invalid laser polarization degree=" << lPolD << G4endl;
    exit(-1);
  }

  useRecoil = true; // switch for electron recoil effect
  if (useRecoil)
  {
    G4cout << "GammaGenerator: electron recoil switched ON!" << G4endl;
  }
  else
  {
    G4cout << "GammaGenerator: electron recoil switched OFF!" << G4endl;
  }
  useNonLin = false; // switch for non-linear interaction effect
  if (useNonLin)
  {
    G4cout << "GammaGenerator: non-linear corrections switched ON!" << G4endl;
  }
  else
  {
    G4cout << "GammaGenerator: non-linear corrections switched OFF!" << G4endl;
  }
  printOut = false; // switch for printing out

  PrintParameters();
}

GammaGenerator::~GammaGenerator()
{
  ComputeIntensity();
  delete fGamma;
}

MCGamma *GammaGenerator::Gamma()
{
  GenerateElectron(); // generate electron in the lab frame
  if (printOut)
    Print4Vector("Initial electron:", Electron4Mom);
  Electron3Boost = Electron4Mom.boostVector();
  GeneratePhoton(); // generate photon in the lab frame
  if (printOut)
    Print4Vector("Initial photon:", Photon4Mom);
  if (useNonLin)
  {
    NonLinEffect(); // apply nonlinear effect to electron 4-momentum
    if (printOut)
      Print4Vector("Dressed electron:", Electron4Mom);
  }
  if (useRecoil && nGamma < 1000)
  { // correction for total cross section
    G4LorentzVector Init4Mom = Electron4Mom + Photon4Mom;
    XSum += (Init4Mom * Init4Mom) / (mec2 * mec2) - 1;
    nGamma++;
  }
  GoToRestFrame(); // move to electron rest frame with aligned photon
  GenerateGamma(); // generate gamma
  if (printOut)
    Print4Vector("At rest gamma:", Gamma4Mom);
  Gamma4Mom = Gamma4Mom.boost(Electron3Boost); // move back to the lab frame
  if (printOut)
    Print4Vector("Final gamma:", Gamma4Mom);
  // now fill the gamma in G4ParticleGun format:
  fGamma->energy = Gamma4Mom.e();
  fGamma->direction = Gamma4Mom.vect().unit();
  fGamma->time = 0. * ns;
  fGamma->position = Electron3Pos; // gamma generated at electron position
  fGamma->polarization = GammaStokes;
  return fGamma;
}

// generate electron momentum 4-vector and position 3-vector in the lab frame
void GammaGenerator::GenerateElectron()
{

  G4double u1 = expRand->fire();
  G4double u2 = expRand->fire();
  G4double u1p = 2 * pi * flatRand->fire();
  G4double u2p = 2 * pi * flatRand->fire();

  G4double x = sqrt(2 * u1 * eEmit * eBeta) * cos(u1p);
  G4double y = sqrt(2 * u2 * eEmit * eBeta) * cos(u2p);
  G4double z = eSigZ * gaussRand->fire(0., 1.);
  G4double xp = -sqrt(2 * u1 * eEmit / eBeta) * (eAlpha * cos(u1p) + sin(u1p));
  G4double yp = -sqrt(2 * u2 * eEmit / eBeta) * (eAlpha * cos(u2p) + sin(u2p));
  G4double eEn = gaussRand->fire(eEner, deEner * eEner) + mec2;
  G4double ePz = sqrt((eEn * eEn - mec2 * mec2) / (1 + xp * xp + yp * yp));
  G4double ePx = ePz * xp;
  G4double ePy = ePz * yp;

  Electron4Mom.setE(eEn);
  Electron4Mom.setPx(ePx);
  Electron4Mom.setPy(ePy);
  Electron4Mom.setPz(ePz);
  Electron3Pos.set(x, y, z);
}

// generate photon momentum 4-vector and polarization 4-vector in the lab frame
void GammaGenerator::GeneratePhoton()
{
  G4double pEn = 2. * pi * hbarc / lLamb;
  G4double elX = Electron3Pos.x();
  G4double elY = Electron3Pos.y();
  G4double elZ = Electron3Pos.z();
  G4double lc1 = elX * elZ / (elZ * elZ + lRayl * lRayl);
  G4double lc2 = elY * elZ / (elZ * elZ + lRayl * lRayl);
  G4double lcX = -lc1 / sqrt(1 + lc1 * lc1 + lc2 * lc2);
  G4double lcY = -lc2 / sqrt(1 + lc1 * lc1 + lc2 * lc2);
  G4double lcZ = -1. / sqrt(1 + lc1 * lc1 + lc2 * lc2);
  Photon4Mom.setE(pEn);
  Photon4Mom.setPx(lcX * pEn);
  Photon4Mom.setPy(lcY * pEn);
  Photon4Mom.setPz(lcZ * pEn);
  Photon4Mom.rotateY(-elCross);
  Photon3Pol.setX(lPolD * cos(lPolA));
  Photon3Pol.setY(lPolD * sin(lPolA));
  Photon3Pol.setZ(0.);
  Photon3Pol.rotateY(-elCross);
}

// apply nonlinear effects to electron 4-momentum and mass
void GammaGenerator::NonLinEffect()
{
  G4double Int = lPowr / (lRayl * lLamb * elFreq * lSigT);
  G4double omega = 2 * pi * c_light / lLamb;
  G4double aL = sqrt(4 * pi * classic_electr_radius * Int * 2 / (omega * omega * mec2 / c_light));
  G4LorentzVector Phot4Versor = Photon4Mom;
  Phot4Versor *= 1. / (Photon4Mom.t());
  G4double corr = 0.5 * pow(aL * mec2, 2) / (Phot4Versor * Electron4Mom);
  Electron4Mom += (Phot4Versor * corr);
  mec2 *= sqrt(1 + aL * aL);
}

// move all initial vectors to electron rest frame and align photon on +z
void GammaGenerator::GoToRestFrame()
{
  G4ThreeVector Electron3Boost0 = Electron4Mom.boostVector();
  Electron3Boost0 *= -1.;
  Electron4Mom = Electron4Mom.boost(Electron3Boost0);
  Photon4Mom = Photon4Mom.boost(Electron3Boost0);
  if (printOut)
    Print4Vector("Electron at rest:", Electron4Mom);
  if (printOut)
    Print4Vector("Photon for electron at rest:", Photon4Mom);
  G4double PhotPhi = Photon4Mom.phi() / rad;
  G4double PhotThe = Photon4Mom.theta() / rad;
  Photon4Mom.rotateZ(-PhotPhi);
  Photon4Mom.rotateY(-PhotThe);
  if (printOut)
    Print4Vector("Photon along z:", Photon4Mom);
}

// generate gamma in electron rest frame
void GammaGenerator::GenerateGamma()
{
  G4double pEn = Photon4Mom.e();
  G4double pPolP = Photon3Pol.phi();
  G4double emin = pEn / (1 + 2 * pEn / mec2);
  G4double emax = pEn;
  G4double gEn = 0. * MeV, gThe = 0. * rad, gPhi = 0. * rad;
  G4double probE = 0., randE = 0., probPhi = 0., randPhi = 0., x = 0., y = 0., Axy = 0., Bxy = 0.;
  do
  {
    gEn = flatRand->fire() * (emax - emin) + emin;
    x = mec2 * (1 / pEn - 1 / gEn);
    y = pEn / gEn + gEn / pEn;
    probE = (x * x + 2 * x + y) / 2 / (1 + pEn / mec2);
    randE = flatRand->fire();
  } while (randE > probE);
  gThe = acos(1 + x);
  do
  {
    gPhi = 2 * pi * flatRand->fire();
    Axy = (x * x + 2 * x + y) / (2 * x * x + 4 * x + y) / pi;
    Bxy = (x * x + 2 * x) / (2 * x * x + 4 * x + y) / pi;
    probPhi = Axy + Bxy * lPolD * cos(2 * pPolP - 2 * gPhi);
    randPhi = flatRand->fire();
  } while (randPhi > probPhi);
  Gamma4Mom.setE(gEn);
  Gamma4Mom.setPx(gEn * sin(gThe) * cos(gPhi));
  Gamma4Mom.setPy(gEn * sin(gThe) * sin(gPhi));
  Gamma4Mom.setPz(gEn * cos(gThe));
  GammaStokes.setX(0.);
  GammaStokes.setY(0.); // circular polarization not implemented
  GammaStokes.setZ(lPolD * (1 - x * x / 2 / (x * x + 2 * x + y)));
}

// compute and print gamma source intensity:
void GammaGenerator::ComputeIntensity()
{
  G4double sCor = 1.;
  if (useRecoil)
    sCor = 1 - XSum / nGamma;
  G4double XSect = 0.665 * sCor * barn;
  G4double eInt = eChar * elFreq;
  G4double scLum = 1. / (2. * pi * eBeta * eEmit + lLamb * lRayl / 2.);
  G4double GammaInt = XSect * scLum * eInt * lPowr * lLamb / (h_Planck * c_light * eplus * elFreq);
  G4cout << G4endl << "Normalize your rates to the gamma source intensity: " << std::scientific << GammaInt / hertz << " gamma/s" << G4endl;
  if (useRecoil)
    G4cout << "Electron recoil correction factor: " << sCor << G4endl;
}

void GammaGenerator::Print4Vector(G4String head, G4LorentzVector v)
{
  G4cout << head << std::scientific << std::setprecision(2) << " E=" << v.e() / MeV << "MeV, Px=" << v.px() / MeV << "MeV, Py=" << v.py() / MeV << "MeV, Pz=" << v.pz() / MeV << "MeV, The=" << v.theta() / mrad << "mrad, Phi=" << v.phi() / deg << "deg" << G4endl;
}

void GammaGenerator::PrintParameters()
{
  G4cout << "GammaGenerator parameters:" << G4endl;
  G4cout << "   Electron beam:" << G4endl;
  G4cout << "   Energy=" << eEner / MeV << "MeV; Energy spread=" << deEner << "; Bunch length=" << eSigZ / mm << "mm; Emittance=" << eEmit / nanometer << "nm; Beta=" << eBeta / m << "m; Number of bunches=" << eNB << "; Bunch charge=" << eChar / coulomb << "C." << G4endl;
  G4cout << "   Laser beam:" << G4endl;
  G4cout << "   Wavelength=" << lLamb / nanometer << "nm; Energy spread=" << dlEner << "; Rayleigh length=" << lRayl / m << "m; Pulse length=" << lSigT / picosecond << "ps; Pulse power=" << lPowr / watt << "W; Polarization degree=" << lPolD << "; Polarization angle=" << lPolA / deg << "deg." << G4endl;
  G4cout << "   Electron-Laser interaction:" << G4endl;
  G4cout << "   Crossing angle=" << elCross / mrad << "mrad; Collision frequency=" << elFreq / megahertz << "MHz." << G4endl << G4endl;
}
