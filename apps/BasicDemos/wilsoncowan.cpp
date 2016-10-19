//Machine-generated by miind.py. Edit at your own risk.

#include <boost/timer/timer.hpp>
#include <GeomLib.hpp>
#include <MPILib/include/MPINetworkCode.hpp>
#include <MPILib/include/RateAlgorithmCode.hpp>
#include <MPILib/include/SimulationRunParameter.hpp>
#include <MPILib/include/report/handler/RootReportHandler.hpp>
#include <MPILib/include/WilsonCowanAlgorithm.hpp>
#include <MPILib/include/PersistantAlgorithm.hpp>
#include <MPILib/include/DelayAlgorithmCode.hpp>
#include <MPILib/include/RateFunctorCode.hpp>

typedef MPILib::MPINetwork<double, MPILib::utilities::CircularDistribution> Network;

int main(int argc, char *argv[]) {
    Network network;
    boost::timer::auto_cpu_timer t;

#ifdef ENABLE_MPI
    // initialise the mpi environment this cannot be forwarded to a class
    boost::mpi::environment env(argc, argv);
#endif

    try {	// generating algorithms
        const MPILib::Time t_mem_0 = 50e-3;
        const double f_noise_0 = 1.0;
        MPILib::Rate f_max_0 = 10.0;
        MPILib::Rate I_ext_0 = 0;
        MPILib::WilsonCowanParameter  par_wil_0(t_mem_0,f_max_0,f_noise_0,I_ext_0);
        MPILib::WilsonCowanAlgorithm alg_wc_0(par_wil_0);
        MPILib::RateAlgorithm<double> rate_alg_1(100.0);
        // generating nodes
        MPILib::NodeId id_0 = network.addNode(alg_wc_0,MPILib::EXCITATORY_GAUSSIAN);
        MPILib::NodeId id_1 = network.addNode(rate_alg_1,MPILib::EXCITATORY_GAUSSIAN);
        // generating connections
        double con_1_0(0.1);
        network.makeFirstInputOfSecond(id_1,id_0,con_1_0);
        // generation simulation parameter
        const MPILib::Time tmin = 0;
        const MPILib::Time tmax = 0.3;
        const MPILib::Rate fmin = 0;
        const MPILib::Rate fmax = 10;
        const MPILib::Potential statemin = 0;
        const MPILib::Potential statemax = 0.02;
        const MPILib::Potential densemin = 0;
        const MPILib::Potential densemax = 250;
        MPILib::CanvasParameter par_canvas(tmin,tmax,fmin,fmax,statemin,statemax,densemin,densemax);

        MPILib::report::handler::RootReportHandler handler("wilsoncowan",true,true, par_canvas);
        handler.addNodeToCanvas(id_0);
        SimulationRunParameter par_run( handler,1000000,0,0.3,1e-03,1e-03,"wilson.log",1e-03);
        network.configureSimulation(par_run);
        network.evolve();
    } catch(std::exception exc) {
        std::cout << exc.what() << std::endl;
#ifdef ENABLE_MPI
        //Abort the MPI environment in the correct way :
        env.abort(1);
#endif
    }

    MPILib::utilities::MPIProxy().barrier();
    t.stop();
    if (MPILib::utilities::MPIProxy().getRank() == 0) {

        std::cout << "Overall time spend\n";
        t.report();
    }
    return 0;
}
