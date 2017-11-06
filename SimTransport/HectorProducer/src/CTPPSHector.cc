#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "SimTransport/HectorProducer/interface/CTPPSHector.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"
#include "HepMC/SimpleVector.h"

#include "CLHEP/Random/RandGauss.h"

#include "TRandom3.h"
#include <TMath.h>
#include <TMatrixT.h>
#include <TH2F.h>

#include "H_Parameters.h"

#include <math.h>

CTPPSHector::CTPPSHector(const edm::ParameterSet & param, bool verbosity,bool CTPPSTransport) : 
    m_smearAng(false),m_sig_e(0.),m_smearE(false),m_sigmaSTX(0.),m_sigmaSTY(0.),m_sigmaSX(0.),m_sigmaSY(0.),
    fCrossAngleCorr(false),fCrossingAngleBeam1(0.),fCrossingAngleBeam2(0.),fBeamMomentum(0),fBeamEnergy(0),
    fVtxMeanX(0.),fVtxMeanY(0.),fVtxMeanZ(0.),fMomentumMin(0.),
    m_verbosity(verbosity), 
    m_CTPPSTransport(CTPPSTransport),NEvent(0)

{
    // Create LHC beam line
    edm::ParameterSet hector_par = param.getParameter<edm::ParameterSet>("CTPPSHector");

    // User definitons
    lengthctpps     = hector_par.getParameter<double>("BeamLineLengthCTPPS" );
    m_f_ctpps_f     = (float) hector_par.getParameter<double>("CTPPSf");
    m_b_ctpps_b     = (float) hector_par.getParameter<double>("CTPPSb");

    beam1filename   = hector_par.getParameter<string>("Beam1");
    beam2filename   = hector_par.getParameter<string>("Beam2");  
    m_smearAng      = hector_par.getParameter<bool>("smearAng");
    m_sigmaSTX      = hector_par.getParameter<double>("sigmaSTX" );
    m_sigmaSTY      = hector_par.getParameter<double>("sigmaSTY" );
    m_sigmaSX       = hector_par.getParameter<double>("sigmaSX");
    m_sigmaSY       = hector_par.getParameter<double>("sigmaSY");
    m_smearE        = hector_par.getParameter<bool>("smearEnergy");
    m_sig_e         = hector_par.getParameter<double>("sigmaEnergy");
    etacut          = hector_par.getParameter<double>("EtaCutForHector" );
    //CTPPS
    fCrossAngleCorr = hector_par.getParameter<bool>("CrossAngleCorr");
    fCrossingAngleBeam1  = hector_par.getParameter<double>("CrossingAngleBeam1");
    fCrossingAngleBeam2  = hector_par.getParameter<double>("CrossingAngleBeam2");
    fBeamEnergy     = hector_par.getParameter<double>("BeamEnergy"); // beam energy in GeV
    fVtxMeanX       = hector_par.getParameter<double>("VtxMeanX");
    fVtxMeanY       = hector_par.getParameter<double>("VtxMeanY");
    fVtxMeanZ       = hector_par.getParameter<double>("VtxMeanZ");
    fMomentumMin    = hector_par.getParameter<double>("MomentumMin"); 

    fBeamMomentum = sqrt(fBeamEnergy*fBeamEnergy - ProtonMassSQ);
    std::vector<std::string> bp = hector_par.getParameter<std::vector<std::string> >("BeamPositions");

    for(std::vector<std::string>::const_iterator it=bp.begin();it!=bp.end();it++) {
       double bz,bx,by;
       std::string::size_type sz,sz2; //idx in *it and substr of *it
       bz=std::stod(*it,&sz);
       bx=std::stod((*it).substr(sz+1),&sz2);
       sz+=sz2+1;
       by=std::stod((*it).substr(sz+1));
       (bz>0)?PosP.push_back(std::make_tuple<double,double,double>((double)bz,(double)bx,(double)by)):
              PosN.push_back(std::make_tuple<double,double,double>((double)(-bz),(double)bx,(double)by));
    }

    theCorrespondenceMap.clear();

    if(m_verbosity) {
        edm::LogInfo("CTPPSHectorSetup") << "===================================================================\n"  
            << " * * * * * * * * * * * * * * * * * * * * * * * * * * * *           \n"  
            << " *                                                         *       \n"  
            << " *                   --<--<--  A fast simulator --<--<--     *     \n"  
            << " *                 | --<--<--     of particle   --<--<--     *     \n"  
            << " *  ----HECTOR----<                                          *     \n"  
            << " *                 | -->-->-- transport through-->-->--      *     \n"   
            << " *                   -->-->-- generic beamlines -->-->--     *     \n"  
            << " *                                                           *     \n"   
            << " * JINST 2:P09005 (2007)                                     *     \n"  
            << " *      X Rouby, J de Favereau, K Piotrzkowski (CP3)         *     \n"  
            << " *       http://www.fynu.ucl.ac.be/hector.html               *     \n"  
            << " *                                                           *     \n"  
            << " * Center for Cosmology, Particle Physics and Phenomenology  *     \n"  
            << " *              Universite catholique de Louvain             *     \n"  
            << " *                 Louvain-la-Neuve, Belgium                 *     \n"  
            << " *                                                         *       \n"  
            << " * * * * * * * * * * * * * * * * * * * * * * * * * * * *           \n"   
            << " CTPPSHector configuration: \n" 
            << " m_CTPPSTransport   = " << m_CTPPSTransport << "\n"
            << " lengthctpps      = " << lengthctpps << "\n"
            << " m_f_ctpps_f      =  " << m_f_ctpps_f << "\n"
            << " m_b_ctpps_b      =  " << m_b_ctpps_b << "\n"
            << "===================================================================\n";
    }  
    edm::FileInPath b1(beam1filename.c_str());
    edm::FileInPath b2(beam2filename.c_str());
    extern int kickers_on;
    kickers_on = 1; // turn on HECTOR kickers if Xangle correction is required

    // construct beam line for CTPPS (forward 1 backward 2):                                                                                           
    if(m_CTPPSTransport && lengthctpps>0. ) {
        m_beamlineCTPPS1 = std::unique_ptr<H_BeamLine>(new H_BeamLine( -1, lengthctpps + 0.1 )); // (direction, length)
        m_beamlineCTPPS2 = std::unique_ptr<H_BeamLine>(new H_BeamLine( 1, lengthctpps + 0.1 )); //
        m_beamlineCTPPS1->fill( b2.fullPath(), 1, "IP5" );
        m_beamlineCTPPS2->fill( b1.fullPath(), 1, "IP5" );
        m_beamlineCTPPS1->offsetElements( 120, 0.097 );
        m_beamlineCTPPS2->offsetElements( 120,-0.097 );
        m_beamlineCTPPS1->calcMatrix();
        m_beamlineCTPPS2->calcMatrix();
        fBeamXatIP=fVtxMeanX*cm_to_mm; // position in mm
        fBeamYatIP=fVtxMeanY*cm_to_mm;
        BeamPositionCalibration(fBeamXatIP,fBeamYatIP,0.);
        BeamProfile();
    } else {
        if ( m_verbosity ) LogDebug("CTPPSHectorSetup") << "CTPPSHector: WARNING: lengthctpps=  " << lengthctpps;
    } 
    kickers_on=0;
}

CTPPSHector::~CTPPSHector(){

    for (std::map<unsigned int,H_BeamParticle*>::iterator it = m_beamPart.begin(); it != m_beamPart.end(); ++it ) {
        delete (*it).second;
    }

}

void CTPPSHector::clearApertureFlags(){
    m_isStoppedctpps.clear();
}

void CTPPSHector::clear(){
    for ( std::map<unsigned int,H_BeamParticle*>::iterator it = m_beamPart.begin(); it != m_beamPart.end(); ++it ) {
        delete (*it).second;
    };
    m_beamPart.clear();
    m_direct.clear();
    m_eta.clear();
    m_pdg.clear();
    m_pz.clear();
    m_isCharged.clear();  
}

void CTPPSHector::add( const HepMC::GenEvent * evt ,const edm::EventSetup & iSetup, CLHEP::HepRandomEngine * engine) {

    H_BeamParticle* h_p  = NULL;
    unsigned int line;

    for (HepMC::GenEvent::particle_const_iterator eventParticle =evt->particles_begin();
            eventParticle != evt->particles_end();
            ++eventParticle ) {
        if ( (*eventParticle)->status() == 1 && (*eventParticle)->pdg_id()==2212 ){
            if ( abs( (*eventParticle)->momentum().eta())>etacut && abs( (*eventParticle)->momentum().pz())>fMomentumMin){
                line = (*eventParticle)->barcode();
                if ( m_beamPart.find(line) == m_beamPart.end() ) {
                    double charge=1.;
                    m_isCharged[line] = false;// neutrals
                    HepMC::GenParticle * g = (*eventParticle);	
                    iSetup.getData( pdt );
                    const ParticleData * part = pdt->particle( g->pdg_id() );
                    if (part){
                        charge = part->charge();
                    }
                    if(charge !=0) m_isCharged[line] = true;//charged
                    double mass = (*eventParticle)->generatedMass();

                    h_p = new H_BeamParticle(mass,charge);

                    double px,py,pz,e;
                    double TXforPosition=0.0, TYforPosition=0.0;//urad

                    px = (*eventParticle)->momentum().px();	  
                    py = (*eventParticle)->momentum().py();	  
                    pz = (*eventParticle)->momentum().pz();	  

                    e = sqrt(pow(mass,2)+pow(px,2)+pow(py,2)+pow(pz,2));

                    TXforPosition=(pz>0)?fCrossingAngleBeam2:fCrossingAngleBeam1; // in the CMS ref. framne
                    // Apply Beam and Crossing Angle Corrections
                    LorentzVector p_out(px,py,pz,e);
                    //LorentzBoost(const_cast<LorentzVector&>(p_out),"LAB");
                    ApplyBeamCorrection(p_out, engine);
                    //TXforPosition=0.;

                    // from mm to cm        
                    double XforPosition = (*eventParticle)->production_vertex()->position().x()/cm;//cm
                    double YforPosition = (*eventParticle)->production_vertex()->position().y()/cm;//cm
                    double ZforPosition = (*eventParticle)->production_vertex()->position().z()/cm;//cm

                    if(m_verbosity) LogDebug("CTPPSHectorEventProcessing") <<
                                    " fVtxMeanX: " << fVtxMeanX << " fVtxMeanY: " << fVtxMeanY << " fVtxMeanZ: "  << fVtxMeanZ ;
                    // It is important to set the Position before the 4Momentum otherwise HECTOR resets variables
                    h_p->setPosition(-((XforPosition-fVtxMeanX)*cm_to_um+fBeamXatIP*mm_to_um),(YforPosition-fVtxMeanY)*cm_to_um+fBeamYatIP*mm_to_um,
                                        -TXforPosition,TYforPosition,-(ZforPosition-fVtxMeanZ)*cm_to_m);
                    
                    std::cout <<"Befor set4momentum -> TX: " << -h_p->getTX() << " " << h_p->getTY() << std::endl;
                    h_p->set4Momentum(-p_out.px(), p_out.py(), abs(p_out.pz()), p_out.e());
                    double tx = atan(-p_out.px()/abs(p_out.pz()))/urad;
                    double tx2 = atan2(-p_out.px(),abs(p_out.pz()))/urad;
                    std::cout <<"After set4momentum -> TX: " << h_p->getTX() << " " << h_p->getTY()
                              << " TX should be:" << -TXforPosition+tx << std::endl;
                    //BeamProfile();
                    m_beamPart[line] = h_p;
                    m_direct[line] = 0;
                    m_direct[line] = ( pz > 0 ) ? 1 : -1;
                    m_eta[line] = p_out.eta();
                    m_pdg[line] = (*eventParticle)->pdg_id();
                    m_pz[line]  = p_out.pz();
                    if(m_verbosity) { 
                        LogDebug("CTPPSHectorEventProcessing") << "CTPPSHector:add: barcode = " << line 
                            << " status = " << g->status() 
                            << " PDG Id = " << g->pdg_id() 
                            << " mass = " << mass 
                            << " pz = " << pz 
                            << " charge = " << charge 
                            << " m_isCharged[line] = " << m_isCharged[line];
                    } 
                }// if find line
            }// if eta > 8.2
        }// if status
    }// for loop

}

void CTPPSHector::filterCTPPS(TRandom3* rootEngine){

    unsigned int line;
    H_BeamParticle * part = NULL;
 
    std::map< unsigned int, H_BeamParticle* >::iterator it;

    bool is_stop;
    int direction;

    float x1_ctpps;
    float y1_ctpps;

    if ( m_beamPart.size() && lengthctpps>0. ) {

        for (it = m_beamPart.begin(); it != m_beamPart.end(); ++it ) {
            line = (*it).first;
            part = (*it).second;

            if(m_verbosity) LogDebug("CTPPSHectorEventProcessing") << "CTPPSHector:filterCTPPS: barcode = " << line;
            if ( (*m_isCharged.find( line )).second ) {
                direction = (*m_direct.find( line )).second;
                if ( direction == 1 && m_beamlineCTPPS1 != 0 ) {
                    
 		    part->computePath(&*m_beamlineCTPPS1);

                    is_stop = part->stopped(&* m_beamlineCTPPS1 );
                    if(m_verbosity) LogDebug("CTPPSHectorEventProcessing") << 
                                    "CTPPSHector:filterCTPPS: barcode = " << line << " positive is_stop=  "<< is_stop;
                }
                else if ( direction == -1 && m_beamlineCTPPS2 != 0 ){

                    part->computePath(&*m_beamlineCTPPS2);

                    is_stop = part->stopped(&*m_beamlineCTPPS2 );
                    if(m_verbosity) LogDebug("CTPPSHectorEventProcessing") << 
                                    "CTPPSHector:filterCTPPS: barcode = " << line << " negative is_stop=  "<< is_stop;
                }
                else {
                    is_stop = true;
                    if(m_verbosity) LogDebug("CTPPSHectorEventProcessing") << "CTPPSHector:filterCTPPS: barcode = " << line << " 0      is_stop=  "<< is_stop;
                }

                //propagating
                m_isStoppedctpps[line] = is_stop;
                if(m_verbosity) LogDebug("CTPPSHectorEventProcessing") << 
                                "CTPPSHector:filterCTPPS: barcode = " << line << " isStopped=" << (*m_isStoppedctpps.find(line)).second;

                if (!is_stop) {
                    if ( direction == 1 ) part->propagate( m_f_ctpps_f ); 
                    if ( direction == -1 ) part->propagate( m_b_ctpps_b );  
                    x1_ctpps = -part->getX()/millimeter;
                    y1_ctpps = part->getY()/millimeter;
                    if(m_verbosity) LogDebug("CTPPSHectorEventProcessing") << 
                                   "CTPPSHector:filterCTPPS: barcode = " << line << " x=  "<< x1_ctpps <<" y= " << y1_ctpps;
                    m_xAtTrPoint[line]  = x1_ctpps;
                    m_yAtTrPoint[line]  = y1_ctpps;
                    m_TxAtTrPoint[line] = -part->getTX(); // needs to be reflected due to the way phi is calculated here
                    m_TyAtTrPoint[line] = part->getTY();
                    m_eAtTrPoint[line]  = part->getE();

                }
            }// if isCharged
            else {
                m_isStoppedctpps[line] = true;// imply that neutral particles stopped to reach 420m
                if(m_verbosity) LogDebug("CTPPSHectorEventProcessing") << 
                                "CTPPSHector:filterCTPPS: barcode = " << line << " isStopped=" << (*m_isStoppedctpps.find(line)).second;
            }

        } // for (it = m_beamPart.begin(); it != m_beamPart.end(); it++ ) 
    } // if ( m_beamPart.size() )

}//

int  CTPPSHector::getDirect( unsigned int part_n ) const {
    std::map<unsigned int, int>::const_iterator it = m_direct.find( part_n );
    if ( it != m_direct.end() ){
        return (*it).second;
    }
    return 0;
}

void CTPPSHector::print() const {
    for (std::map<unsigned int,H_BeamParticle*>::const_iterator it = m_beamPart.begin(); it != m_beamPart.end(); ++it ) {
        (*it).second->printProperties();
    };
}

void CTPPSHector::ApplyBeamCorrection(LorentzVector& p_out, CLHEP::HepRandomEngine* engine)
{

    double theta = p_out.theta();
    double thetax = atan(p_out.px()/p_out.pz());
    double thetay = atan(p_out.py()/p_out.pz());
    double energy = p_out.e();

    int direction = (p_out.pz()>0)?1:-1;

    if (p_out.pz()<0) theta=CLHEP::pi-theta;

    double dtheta_x = (double)(m_smearAng)?CLHEP::RandGauss::shoot(engine,0.,m_sigmaSTX):0;
    double dtheta_y = (double)(m_smearAng)?CLHEP::RandGauss::shoot(engine,0.,m_sigmaSTY):0;
    double denergy  = (double)(m_smearE)?CLHEP::RandGauss::shoot(engine,0.,m_sig_e):0.;

    double s_theta = sqrt(pow(thetax+dtheta_x*urad,2)+pow(thetay+dtheta_y*urad,2)); 
    double s_phi = atan2(thetay+dtheta_y*urad,thetax+dtheta_x*urad);
    energy+=denergy;
    double p = sqrt(pow(energy,2)-ProtonMassSQ);

    p_out.setPx((double)p*sin(s_theta)*cos(s_phi));
    p_out.setPy((double)p*sin(s_theta)*sin(s_phi));
    p_out.setPz((double)p*(cos(s_theta))*direction);
    p_out.setE(energy);
}

void CTPPSHector::LorentzBoost(H_BeamParticle& h_p, const string& frame)
{
     LorentzVector p_out = HectorParticle2LorentzVector(h_p);
     LorentzBoost(p_out,frame);
     h_p.set4Momentum(-p_out.px(),p_out.py(),p_out.pz(),p_out.e());
}
void CTPPSHector::LorentzBoost(LorentzVector& p_out, const string& frame)
{
    // Use a matrix
    double microrad = 1.e-6;
    TMatrixT<Double_t> tmpboost(4,4);
    double alpha_ = 0.;
    double fBoostAngle = (p_out.pz()>0)?fCrossingAngleBeam2*microrad:-fCrossingAngleBeam1*microrad;
    tmpboost(0,0) = 1./cos(fBoostAngle);
    tmpboost(0,1) = -cos(alpha_)*sin(fBoostAngle);
    tmpboost(0,2) = -tan(fBoostAngle)*sin(fBoostAngle);
    tmpboost(0,3) = -sin(alpha_)*sin(fBoostAngle);
    tmpboost(1,0) = -cos(alpha_)*tan(fBoostAngle);
    tmpboost(1,1) = 1.;
    tmpboost(1,2) = cos(alpha_)*tan(fBoostAngle);
    tmpboost(1,3) = 0.;
    tmpboost(2,0) = 0.;
    tmpboost(2,1) = -cos(alpha_)*sin(fBoostAngle);
    tmpboost(2,2) = cos(fBoostAngle);
    tmpboost(2,3) = -sin(alpha_)*sin(fBoostAngle);
    tmpboost(3,0) = -sin(alpha_)*tan(fBoostAngle);
    tmpboost(3,1) = 0.;
    tmpboost(3,2) = sin(alpha_)*tan(fBoostAngle);
    tmpboost(3,3) = 1.;

    if(frame=="LAB") tmpboost.Invert();

    TMatrixT<Double_t> p4(4,1);
    p4(0,0) = p_out.e();
    p4(1,0) = p_out.px();
    p4(2,0) = p_out.py();
    p4(3,0) = p_out.pz();
    TMatrixT<Double_t> p4lab(4,1);
    p4lab = tmpboost * p4;
    p_out.setPx(p4lab(1,0));
    p_out.setPy(p4lab(2,0));
    p_out.setPz(p4lab(3,0));
    p_out.setE(p4lab(0,0));
}

HepMC::GenEvent * CTPPSHector::addPartToHepMC( HepMC::GenEvent * evt ){
    NEvent++;
    theCorrespondenceMap.clear();

    unsigned int line;

    HepMC::GenParticle * gpart;
    long double tx,ty,theta,fi,energy,time = 0;
    std::map< unsigned int, H_BeamParticle* >::iterator it;

    for (it = m_beamPart.begin(); it != m_beamPart.end(); ++it ) {
        line = (*it).first;
        if(!m_CTPPSTransport) m_isStoppedctpps[line] = true;
        if(m_verbosity) {
            LogDebug("CTPPSHectorEventProcessing") << "CTPPSHector:addPartToHepMC: barcode = " << line << "\n"
                << "CTPPSHector:addPartToHepMC: isStoppedctpps=" << (*m_isStoppedctpps.find(line)).second;
        }
        if (!((*m_isStoppedctpps.find(line)).second)){

            gpart = evt->barcode_to_particle( line );
            if ( gpart ) {
                tx     = (*m_TxAtTrPoint.find(line)).second / 1000000.;
                ty     = (*m_TyAtTrPoint.find(line)).second / 1000000.;
                theta  = sqrt((tx*tx) + (ty*ty));
                double ddd = 0.;
                long double fi_  = 0.; 
                if ( !((*m_isStoppedctpps.find(line)).second)) {
                    if( (*m_direct.find( line )).second >0 ) {
                        ddd = m_f_ctpps_f;
                        fi_    = std::atan2(tx,ty); // tx, ty never == 0?
                    }
                    else if((*m_direct.find( line )).second <0 ) {
                        ddd = m_b_ctpps_b;
                        theta= CLHEP::pi-theta;
                        fi_    = std::atan2(tx,ty); // tx, ty never == 0?
                    }
                } 
                fi = fi_; 
                energy = (*m_eAtTrPoint.find(line)).second;
                time = ( ddd*meter - gpart->production_vertex()->position().z()*mm ); // mm

                if(ddd != 0.) {
                    if(m_verbosity) {
                        LogDebug("CTPPSHectorEventProcessing") <<"CTPPSHector:: x= "<< (*(m_xAtTrPoint.find(line))).second*0.001<< "\n"
                            <<"CTPPSHector:: y= "<< (*(m_yAtTrPoint.find(line))).second*0.001<< "\n"
                            <<"CTPPSHector:: z= "<< ddd * (*(m_direct.find( line ))).second*1000. << "\n"
                            <<"CTPPSHector:: t= "<< time;
                    }

                    HepMC::GenVertex * vert = new HepMC::GenVertex( HepMC::FourVector( (*(m_xAtTrPoint.find(line))).second*0.001,
                                (*(m_yAtTrPoint.find(line))).second*0.001,
                                ddd * (*(m_direct.find( line ))).second*1000.,
                                time + .001*time ) );

                    gpart->set_status( 2 );
                    vert->add_particle_in( gpart );
                    double Pmom = sqrt(energy*energy - ProtonMassSQ);
                    vert->add_particle_out( new HepMC::GenParticle( HepMC::FourVector(Pmom*std::sin(theta)*std::sin(fi),
                                    Pmom*std::sin(theta)*std::cos(fi),
                                    Pmom*std::cos(theta),
                                    energy ),gpart->pdg_id(), 1, gpart->flow() ) );
                    evt->add_vertex( vert );

                    int ingoing = (*vert->particles_in_const_begin())->barcode();
                    int outgoing = (*vert->particles_out_const_begin())->barcode();
                    LHCTransportLink theLink(ingoing,outgoing);
                    if (m_verbosity) LogDebug("CTPPSHectorEventProcessing") << "CTPPSHector:addPartToHepMC: LHCTransportLink " << theLink;
                    theCorrespondenceMap.push_back(theLink);

                    if(m_verbosity) LogDebug("CTPPSHectorEventProcessing") << "CTPPSHector::TRANSPORTED pz= " << gpart->momentum().pz()  
                        << " eta= "<< gpart->momentum().eta() << " status= "<< gpart->status();
                }// ddd
            }// if gpart
        }// if !isStopped

        else {
            gpart = evt->barcode_to_particle( line );
            if ( gpart ) {
                gpart->set_status( 2 );
                if(m_verbosity) LogDebug("CTPPSHectorEventProcessing") << "CTPPSHector::NON-transp. pz= " << gpart->momentum().pz()  
                    << " eta= "<< gpart->momentum().eta()  << " status= "<< gpart->status();
            }
        }
    }//for 

    return evt;
} 
void CTPPSHector::BeamPositionCalibration(double& xpos,double& ypos,double q2)
{
     double sigx_target=1e-3;
     double sigy_target=1e-3;
//
     double deltax=(q2>0)?0.0:1.0;
     double deltay=1.0;
     double sigx = 0.;
     double sigy = 0.;
     double sigx_min = 99.;
     double sigy_min = 99.;
     double last_sigx=0.;
     double last_sigy=0.;
     H_BeamParticle h_pp;
     H_BeamParticle h_pn;
//
     //LorentzBoost(h_pp,"LAB");
     //LorentzBoost(h_pn,"LAB");
//
     int dirx=1;
     int diry=1;
     int nInteractions=0;
     int nDirchgx = 0; // number of direction changing with current delta
     int nDirchgy = 0; // number of direction changing with current delta
     while(1) {
        nInteractions++;
        h_pp.setPosition(-xpos*mm_to_um,ypos*mm_to_um,h_pp.getTX()-fCrossingAngleBeam2,0,-fVtxMeanZ*cm_to_m); // the position is given in the CMS frame
        h_pn.setPosition(-xpos*mm_to_um,ypos*mm_to_um,h_pn.getTX()-fCrossingAngleBeam1,0,-fVtxMeanZ*cm_to_m);
        //h_pp.setPosition(-xpos*mm_to_um,ypos*mm_to_um,h_pp.getTX(),0,-fVtxMeanZ*cm_to_m); // the position is given in the CMS frame
        //h_pn.setPosition(-xpos*mm_to_um,ypos*mm_to_um,h_pn.getTX(),0,-fVtxMeanZ*cm_to_m);
        h_pp.computePath(&*m_beamlineCTPPS1);
        h_pn.computePath(&*m_beamlineCTPPS2);
        sigx=0.;
        sigy=0.;
        BdistP.clear();
        BdistN.clear();
        for(unsigned int i=0;i<PosP.size();i++) {
           h_pp.propagate(std::get<0>(PosP.at(i)));
           h_pn.propagate(std::get<0>(PosN.at(i)));

           double xp = std::get<1>(PosP.at(i)); double yp = std::get<2>(PosP.at(i)); // the beam position is given as positive (LHC frame)
           double xn = std::get<1>(PosN.at(i)); double yn = std::get<2>(PosN.at(i)); // ibdem

           BdistP.push_back(make_tuple<double,double,double>((double)std::get<0>(PosP.at(i)),h_pp.getX()*um_to_mm,h_pp.getY()*um_to_mm));
           BdistN.push_back(make_tuple<double,double,double>((double)std::get<0>(PosN.at(i)),h_pn.getX()*um_to_mm,h_pn.getY()*um_to_mm));
       
           sigx+=(pow((h_pp.getX()*um_to_mm-xp),2)+pow((h_pn.getX()*um_to_mm-xn),2));
           sigy+=(pow((h_pp.getY()*um_to_mm-yp),2)+pow((h_pn.getY()*um_to_mm-yn),2));
        }
        sigx=sqrt(sigx);
        sigy=sqrt(sigy);
        if (sigx<sigx_min) {
           sigx_min=sigx; // good, go on in this path
        } else {
           if (sigx>last_sigx) { dirx*=-1;nDirchgx++;}  // change direction
           if (nInteractions>1) {deltax*=0.9;nDirchgx=0;} // decrease delta
        }
        if (sigy<sigy_min) {
           sigy_min=sigy;
        } else {
           if (sigy>last_sigy) {diry*=-1;nDirchgy++;}
           if (nInteractions>1) {deltay*=0.9;nDirchgy=0;}
        }
        double rel_diff_x =(q2>0)?0.:sigx_target-sigx;
        double rel_diff_y =sigy_target-sigy;
        if (rel_diff_x>0&&rel_diff_y>0) break;
        if (deltax==0&&deltay==0) break;
        last_sigx=sigx;
        last_sigy=sigy;
        xpos-=(dirx*deltax);
        ypos-=(diry*deltay);
     }
     //if (m_verbosity){
        //LogDebug("CTPPSHector::BeamPositionCalibration") 
          std::cout
               << "Interaction number " << nInteractions << "\tX = " << xpos << " (mm) \tSigmaX = " << sigx_min << "\tDeltaX = " << deltax << "\n"
               << "                   "                  << "\tY = " << ypos << " (mm) \tSigmaY = " << sigy_min << "\tDeltaY = " << deltay << "\n"
               << "Calibrated beam positions:" << "\n"
               << " Z (m)   \t X (twiss) \t X (calib) \t Y (twiss) \t Y (calib) \t Delta X \t Delta Y\n";
        for(unsigned int i=0;i<BdistP.size();i++) {
            //LogDebug("CTPPSHector::BeamPositionCalibration")
                std::cout 
                <<  std::get<0>(BdistP.at(i))<< " \t "<< std::get<1>(PosP.at(i))<< " \t "<< std::get<1>(BdistP.at(i))
                                              << " \t "<< std::get<2>(PosP.at(i))<< " \t " << std::get<2>(BdistP.at(i))
                                              << " \t "<< std::get<1>(BdistP.at(i))-std::get<1>(PosP.at(i))
                                              << " \t "<< std::get<2>(BdistP.at(i))-std::get<2>(PosP.at(i))
                                              << "\n";
        }
        for(unsigned int i=0;i<BdistN.size();i++) {
           //LogDebug("CTPPSHector::BeamPositionCalibration")
                std::cout 
                 << -std::get<0>(BdistN.at(i))<< " \t "<<std::get<1>(PosN.at(i)) << " \t "<< std::get<1>(BdistN.at(i))
                                              << " \t "<<std::get<2>(PosN.at(i)) << " \t "<<std::get<2>(BdistN.at(i))
                                              << " \t "<<std::get<1>(BdistN.at(i))-std::get<1>(PosN.at(i))
                                              << " \t "<<std::get<2>(BdistN.at(i))-std::get<2>(PosN.at(i))
                                              << "\n";
        }
     //}
     return;
}
void CTPPSHector::BeamProfile()
{    
     TH2F* bp1f = new TH2F("bp1f","",1000,-10000,10000,1000,-10000,10000);
     TH2F* bp2f = new TH2F("bp2f","",1000,-10000,10000,1000,-10000,10000);
     TH2F* bp3f = new TH2F("bp3f","",1000,-10000,10000,1000,-10000,10000);
     TH2F* bp1b = new TH2F("bp1b","",1000,-10000,10000,1000,-10000,10000);
     TH2F* bp2b = new TH2F("bp2b","",1000,-10000,10000,1000,-10000,10000);
     TH2F* bp3b = new TH2F("bp3b","",1000,-10000,10000,1000,-10000,10000);
     for(int j=0;j<2;j++) {
        std::vector<std::tuple<double,double,double> >* DetPos;
        TH2F* prof1;TH2F* prof2; TH2F* prof3;
        H_BeamLine* beamline;
        double crang=0.;
        switch(j) {
              case 0:  // positive side
                     prof1 = &(*bp1f);prof2 = &(*bp2f);prof3 = &(*bp3f);
                     beamline= &(*m_beamlineCTPPS1);   // CTPPS1 corresponds to beam2
                     DetPos = &(PosP);
                     crang=fCrossingAngleBeam2;
                     break;
              case 1: // negative side
                     prof1 = &(*bp1b);prof2 = &(*bp2b); prof3 = &(*bp3b);
                     beamline= &(*m_beamlineCTPPS2);   // CTPPS2 corresponds to beam1
                     DetPos = &(PosN);
                     crang=+fCrossingAngleBeam1;
                     break;
        }
        //crang=0.;
        for(int i=0;i<2000;i++) {
           H_BeamParticle h_p; // Hector always gives a positive pz
           h_p.setPosition(fBeamXatIP*mm_to_um,fBeamYatIP*mm_to_um,h_p.getTX()-crang,h_p.getTY(),-fVtxMeanZ*cm_to_m);
           h_p.smearPos(m_sigmaSX,m_sigmaSY); h_p.smearAng(m_sigmaSTX,m_sigmaSTY); h_p.smearE(m_sig_e);

           //LorentzBoost(h_p,"LAB");

           h_p.computePath(beamline);
           h_p.propagate(std::get<0>(DetPos->at(0)));
           prof1->Fill(h_p.getX(),h_p.getY());
           h_p.propagate(std::get<0>(DetPos->at(1)));
           prof2->Fill(h_p.getX(),h_p.getY());
           h_p.propagate(std::get<0>(DetPos->at(2)));
           prof3->Fill(h_p.getX(),h_p.getY());
        }
     }
     double _beamX_Det1_f    = bp1f->GetMean(1); double _beamY_Det1_f    = bp1f->GetMean(2);
     double _beamSigX_Det1_f = bp1f->GetRMS(1);  double _beamSigY_Det1_f = bp1f->GetRMS(2);
     double _beamX_Det2_f    = bp2f->GetMean(1); double _beamY_Det2_f    = bp2f->GetMean(2);
     double _beamSigX_Det2_f = bp2f->GetRMS(1);  double _beamSigY_Det2_f = bp2f->GetRMS(2);
     double _beamX_Det3_f    = bp3f->GetMean(1); double _beamY_Det3_f    = bp3f->GetMean(2);
     double _beamSigX_Det3_f = bp3f->GetRMS(1);  double _beamSigY_Det3_f = bp3f->GetRMS(2);
     double _beamX_Det1_b    = bp1b->GetMean(1); double _beamY_Det1_b    = bp1b->GetMean(2);
     double _beamSigX_Det1_b = bp1b->GetRMS(1);  double _beamSigY_Det1_b = bp1b->GetRMS(2);
     double _beamX_Det2_b    = bp2b->GetMean(1); double _beamY_Det2_b    = bp2b->GetMean(2);
     double _beamSigX_Det2_b = bp2b->GetRMS(1);  double _beamSigY_Det2_b = bp2b->GetRMS(2);
     double _beamX_Det3_b    = bp3b->GetMean(1); double _beamY_Det3_b    = bp3b->GetMean(2);
     double _beamSigX_Det3_b = bp3b->GetRMS(1);  double _beamSigY_Det3_b = bp3b->GetRMS(2);
     //if (m_verbosity){
         //LogDebug("CTPPSHector::BeamProfile")
           std::cout
             << "HectorForCTPPS: BEAM parameters:\n"
             << "Beam position at Det1 positive side --> " << _beamX_Det1_f << "("<< _beamSigX_Det1_f<<"),\t"<< _beamY_Det1_f << "("<< _beamSigY_Det1_f<<")\n"
             << "Beam position at Det2 positive size --> " << _beamX_Det2_f << "("<< _beamSigX_Det2_f<<"),\t"<< _beamY_Det2_f << "("<< _beamSigY_Det2_f<<")\n"
             << "Beam position at ToF  positive size --> " << _beamX_Det3_f << "("<< _beamSigX_Det3_f<<"),\t"<< _beamY_Det3_f << "("<< _beamSigY_Det3_f<<")\n"
             << "Beam position at Det1 negative side --> " << _beamX_Det1_b << "("<< _beamSigX_Det1_b<<"),\t"<< _beamY_Det1_b << "("<< _beamSigY_Det1_b<<")\n"
             << "Beam position at Det2 negative size --> " << _beamX_Det2_b << "("<< _beamSigX_Det2_b<<"),\t"<< _beamY_Det2_b << "("<< _beamSigY_Det2_b<<")\n"
             << "Beam position at ToF  negative size --> " << _beamX_Det3_b << "("<< _beamSigX_Det3_b<<"),\t"<< _beamY_Det3_b << "("<< _beamSigY_Det3_b<<")\n"
             << "\nBeam positions displacement to the closed orbit:\n"
             << "at Det1 positive side  --> X= " << _beamX_Det1_f*um_to_mm-std::get<1>(PosP.at(0)) << "\tY= "<< _beamY_Det1_f*um_to_mm-std::get<2>(PosP.at(0))<< "\n"
             << "at Det2 positive side  --> X= " << _beamX_Det2_f*um_to_mm-std::get<1>(PosP.at(1)) << "\tY= "<<_beamY_Det2_f*um_to_mm-std::get<2>(PosP.at(1))<< "\n"
             << "at ToF  positive side  --> X= " << _beamX_Det3_f*um_to_mm-std::get<1>(PosP.at(2)) << "\tY= "<<_beamY_Det3_f*um_to_mm-std::get<2>(PosP.at(2))<<"\n"
             << "at Det1 negative side  --> X= " << _beamX_Det1_b*um_to_mm-std::get<1>(PosN.at(0)) << "\tY= "<< _beamY_Det1_b*um_to_mm-std::get<2>(PosN.at(0))<<"\n"
             << "at Det2 negative side  --> X= " << _beamX_Det2_b*um_to_mm-std::get<1>(PosN.at(1)) << "\tY= "<< _beamY_Det2_b*um_to_mm-std::get<2>(PosN.at(1))<<"\n"
             << "at ToF  negative side  --> X= " << _beamX_Det3_b*um_to_mm-std::get<1>(PosN.at(2)) << "\tY= "<< _beamY_Det3_b*um_to_mm-std::get<2>(PosN.at(2))<<"\n";
     //} 
     delete bp1f;delete bp2f;delete bp3f;
             delete bp1b;delete bp2b;delete bp3b;
}
CLHEP::HepLorentzVector CTPPSHector::HectorParticle2LorentzVector(H_BeamParticle hp)
{
     double partP = sqrt(pow(hp.getE(),2)-ProtonMassSQ);
     double theta = sqrt(pow(hp.getTX(),2)+pow(hp.getTY(),2))*urad;
     double pz = partP*cos(theta);
     double px = -tan((double)hp.getTX()*urad)*pz;//PartP*sin(theta)*cos(phi);
     double py = tan((double)hp.getTY()*urad)*pz;//partP*sin(theta)*sin(phi);
     return LorentzVector(px,py,pz,hp.getE());
}
