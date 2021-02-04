#include "MiindOpenSim.h"
#include "NeuralController.h"

using namespace SimTK;
using namespace OpenSim;

MiindOpenSim::MiindOpenSim()
: table_reporter(new TableReporter()),
  reporter(new ConsoleReporter()),
  neural_controller(new NeuralController("neural.mmxml")) {
    // Register Afferent Muscle Type so it can be read by the XML
    Object::registerType(Millard12EqMuscleWithAfferents());
    //model = new Model("MoBL_ARMS_module2_4_onemuscle.osim");
    model = new Model("MoBL_ARMS_module2_4_onemuscle_afferent.osim");
}

MiindOpenSim::~MiindOpenSim() {
    if (table_reporter)
        delete table_reporter;
    if (reporter)
        delete reporter;
    if (integrator)
        delete integrator;
    if (time_stepper)
        delete time_stepper;
    if (neural_controller)
        delete neural_controller;
    if (model)
        delete model;
}

void MiindOpenSim::setElbowExtensionPosture() {
    model->updCoordinateSet().get("elv_angle").setDefaultLocked(true);
    model->updCoordinateSet().get("elv_angle").setDefaultValue(0.0);

    model->updCoordinateSet().get("shoulder_elv").setDefaultLocked(true);
    model->updCoordinateSet().get("shoulder_elv").setDefaultValue(0.0);

    model->updCoordinateSet().get("shoulder_rot").setDefaultLocked(true);
    model->updCoordinateSet().get("shoulder_rot").setDefaultValue(0.0);

    model->updCoordinateSet().get("elbow_flexion").setDefaultLocked(false);
    model->updCoordinateSet().get("elbow_flexion").setDefaultValue(0.298132); // 0.698132 = 40 degrees in radians

    model->updCoordinateSet().get("pro_sup").setDefaultLocked(true);
    model->updCoordinateSet().get("pro_sup").setDefaultValue(0.0);

    model->updCoordinateSet().get("deviation").setDefaultLocked(true);
    model->updCoordinateSet().get("deviation").setDefaultValue(0.0);

    model->updCoordinateSet().get("flexion").setDefaultLocked(true);
    model->updCoordinateSet().get("flexion").setDefaultValue(0.0);
}

void MiindOpenSim::buildSimulation() {
    model->setUseVisualizer(true);

    setElbowExtensionPosture();

    auto biceps_long = (Millard12EqMuscleWithAfferents*)&model->updMuscles().get("BIClong");
    biceps_long->setupAfferents();

    auto triceps_long = (OpenSim::Millard12EqMuscleWithAfferents*)&model->updMuscles().get("TRIlong");
    triceps_long->setupAfferents();

    auto triceps_med = (OpenSim::Millard12EqMuscleWithAfferents*)&model->updMuscles().get("TRImed");
    triceps_med->setupAfferents();

    auto delt_ant = (OpenSim::Millard12EqMuscleWithAfferents*)&model->updMuscles().get("DELT1");
    delt_ant->setupAfferents();

    auto delt_med = (OpenSim::Millard12EqMuscleWithAfferents*)&model->updMuscles().get("DELT2");
    delt_med->setupAfferents();

    auto delt_post = (OpenSim::Millard12EqMuscleWithAfferents*)&model->updMuscles().get("DELT3");
    delt_post->setupAfferents();

    auto fcr = (OpenSim::Millard12EqMuscleWithAfferents*)&model->updMuscles().get("FCR");
    fcr->setupAfferents();

    auto ecr_l = (OpenSim::Millard12EqMuscleWithAfferents*)&model->updMuscles().get("ECRL");
    ecr_l->setupAfferents();

    auto pec_m_1 = (OpenSim::Millard12EqMuscleWithAfferents*)&model->updMuscles().get("PECM1");
    pec_m_1->setupAfferents();

    auto pec_m_2 = (OpenSim::Millard12EqMuscleWithAfferents*)&model->updMuscles().get("PECM2");
    pec_m_2->setupAfferents();

    auto pec_m_3 = (OpenSim::Millard12EqMuscleWithAfferents*)&model->updMuscles().get("PECM3");
    pec_m_3->setupAfferents();

    auto lat_1 = (OpenSim::Millard12EqMuscleWithAfferents*)&model->updMuscles().get("LAT1");
    lat_1->setupAfferents();

    auto lat_2 = (OpenSim::Millard12EqMuscleWithAfferents*)&model->updMuscles().get("LAT2");
    lat_2->setupAfferents();

    auto lat_3 = (OpenSim::Millard12EqMuscleWithAfferents*)&model->updMuscles().get("LAT3");
    lat_3->setupAfferents();

    reporter->set_report_time_interval(0.1);
    reporter->addToReport(biceps_long->getSpindle()->getOutput("primary_Ia"));
    reporter->addToReport(biceps_long->getSpindle()->getOutput("secondary_II"));
    reporter->addToReport(biceps_long->getGTO()->getOutput("gto_out"));
    reporter->addToReport(biceps_long->getOutput("fiber_force"));
    reporter->addToReport(biceps_long->getOutput("fiber_length"));
    reporter->addToReport(biceps_long->getOutput("fiber_velocity"));
    reporter->addToReport(biceps_long->getOutput("lpf_velocity"));
    reporter->addToReport(biceps_long->getOutput("lpf_acceleration"));
    model->addComponent(reporter);

    table_reporter->set_report_time_interval(0.01);
    table_reporter->addToReport(biceps_long->getSpindle()->getOutput("primary_Ia"));
    table_reporter->addToReport(biceps_long->getSpindle()->getOutput("secondary_II"));
    table_reporter->addToReport(biceps_long->getGTO()->getOutput("gto_out"));
    table_reporter->addToReport(biceps_long->getOutput("fiber_force"));
    table_reporter->addToReport(biceps_long->getOutput("fiber_length"));
    table_reporter->addToReport(biceps_long->getOutput("lpf_velocity"));
    table_reporter->addToReport(biceps_long->getOutput("lpf_acceleration"));
    model->addComponent(table_reporter);

    neural_controller->setActuators(model->updActuators());

    neural_controller->prescribeControlForActuator("BIClong", "BIClong_alpha", "BIClong_beta", "BIClong_gamma");
    neural_controller->prescribeControlForActuator("TRIlong", "TRIlong_alpha", "TRIlong_beta", "TRIlong_gamma");
    neural_controller->prescribeControlForActuator("TRImed", "TRImed_alpha", "TRImed_beta", "TRImed_gamma");
    neural_controller->prescribeControlForActuator("DELT1", "DELT1_alpha", "DELT1_beta", "DELT1_gamma");
    neural_controller->prescribeControlForActuator("DELT2", "DELT2_alpha", "DELT2_beta", "DELT2_gamma");
    neural_controller->prescribeControlForActuator("DELT3", "DELT3_alpha", "DELT3_beta", "DELT3_gamma");
    neural_controller->prescribeControlForActuator("FCR", "FCR_alpha", "FCR_beta", "FCR_gamma");
    neural_controller->prescribeControlForActuator("ECRL", "ECRL_alpha", "ECRL_beta", "ECRL_gamma");
    neural_controller->prescribeControlForActuator("PECM1", "PECM1_alpha", "PECM1_beta", "PECM1_gamma");
    neural_controller->prescribeControlForActuator("PECM2", "PECM2_alpha", "PECM2_beta", "PECM2_gamma");
    neural_controller->prescribeControlForActuator("PECM3", "PECM3_alpha", "PECM3_beta", "PECM3_gamma");
    neural_controller->prescribeControlForActuator("LAT1", "LAT1_alpha", "LAT1_beta", "LAT1_gamma");
    neural_controller->prescribeControlForActuator("LAT2", "LAT2_alpha", "LAT2_beta", "LAT2_gamma");
    neural_controller->prescribeControlForActuator("LAT3", "LAT3_alpha", "LAT3_beta", "LAT3_gamma");
    model->addController(neural_controller);
}

void MiindOpenSim::beginSimulation() {
    SimTK::State& state = model->initSystem();

    // Configure the visualizer.
    SimTK::Visualizer& viz = model->updVisualizer().updSimbodyVisualizer();
    viz.setBackgroundType(viz.SolidColor);
    viz.setBackgroundColor(SimTK::Black);
    viz.setShutdownWhenDestructed(true);

    model->equilibrateMuscles(state);

    // begin simulation
    StatesTrajectory states;
    integrator = new SimTK::RungeKutta3Integrator(model->getSystem());
    integrator->setFixedStepSize(0.001);
    integrator->setAccuracy(1e-006);
    integrator->setAllowInterpolation(true);
    integrator->setProjectInterpolatedStates(true);
    time_stepper = new SimTK::TimeStepper(model->getSystem(), *integrator);
    time_stepper->initialize(state);
    time_stepper->setReportAllSignificantStates(true);
    integrator->setReturnEveryInternalStep(true);
    std::cout << "Starting Simulation...\n";
}

void MiindOpenSim::stepSimulation() {
    StatesTrajectory states;
    std::cout << "(" << time_stepper->getState().getTime() << ") " << integrator->getSuccessfulStepStatusString(time_stepper->stepTo(neural_controller->getSimulationTime())) << "                  \r" << std::flush;
    states.append(time_stepper->getState());
}

void MiindOpenSim::finishSimulation() {
    model->updVisualizer().updSimbodyVisualizer().shutdown();

    auto table = table_reporter->getTable();

    ofstream dataFile;
    dataFile.open("output.txt");
    for (int i = 0; i < table.getNumRows(); i++) {
        for (int j = 0; j < 7; j++) { // Currently recording 7 metrics...
            dataFile << table.getRowAtIndex(i).getAsVector()[j] << "," << std::flush;
        }
        dataFile << "\n" << std::flush;
    }
    dataFile.close();
}