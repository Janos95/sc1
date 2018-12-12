#include "assemble_problem.hpp"
#include "constant_size_matrix.hpp"

#include <iostream>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>

#include <unsupported/Eigen/IterativeSolvers>

#include <fmt/core.h>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <boost/timer/timer.hpp>

#include <pcl/visualization/pcl_visualizer.h>


pcl::PolygonMesh toMesh(const Eigen::VectorXd& u, int n)
{
    pcl::PointCloud<pcl::PointXYZ> cloud(u.size(), 1);
    std::vector<::pcl::Vertices > polygons;

    const auto mapToIndex = [=](const int x, const int y){ return x  + y * (n - 1);};
    float meshSize = n;

    for(int x = 0; x < n-1; ++x)
    {
        for(int y = 0; y < 2*n-1; ++y)
        {
            int i = mapToIndex(x,y);
            assert(i < u.size());
            pcl::PointXYZ p((x+1)/meshSize, (y+1)/meshSize, (float)u[i]);
            cloud[i] = p;

            if(x < n-2 && y < 2*n-2)
            {

                pcl::Vertices triangle1, triangle2;
                triangle1.vertices = {uint32_t(i),uint32_t(i+1),uint32_t(i+n-1+1)};
                triangle2.vertices = {uint32_t(i),uint32_t(i+n-1),uint32_t(i+n-1+1)};

                polygons.push_back(std::move(triangle1));
                polygons.push_back(std::move(triangle2));
            }

        }
    }

    pcl::PolygonMesh mesh;
    mesh.polygons = polygons;
    pcl::toPCLPointCloud2(cloud,mesh.cloud);

    return mesh;

}


namespace po = boost::program_options;
namespace fs = boost::filesystem;

int main(int argc, char** argv) {
    boost::timer::auto_cpu_timer all("Whole pipeline took %w seconds\n");

    int N = 100;
    int maxIterations = 10000;
    double tolerance = 1e-10;

    po::options_description desc("Allowed options");
    desc.add_options()
            ("help,h", "produce help message")
            ("mesh-size", po::value<int>()->required(), "size of discretisation")
            ("max-iterations", po::value<int>(), "maximum number of iterations in cg")
            ("tolerance", po::value<double>(), "tolerance in cg");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    bool badargs = false;
    try { po::notify(vm); }
    catch (...) { badargs = true; }

    if (vm.count("help") || badargs) {
        fmt::print("\n\nThis program computes a solution to -laplace(u)=1 and u = 0 on the boundary.\n");
        fmt::print("The domain is (0,2)x(0,1)\n\n");
        fmt::print("Example usage: {} --mesh-size 500 [OPTS]\n\n", fs::basename(argv[0]));
        //fmt::print(std::cout, "{}\n", desc);
        std::cout << desc << std::endl;
        return 1;
    }

    if (vm.count("mesh-size"))
        N = vm["mesh-size"].as<int>();
    if (vm.count("tolerance"))
        tolerance = vm["tolerance"].as<double>();
    if (vm.count("max-iterations"))
        maxIterations = vm["max-iterations"].as<int>();


    Eigen::VectorXd b = Eigen::VectorXd::Constant((2 * N - 1) * (N - 1), 1.0) * 1.0 / (N * N);
    ConstantSizeMatrix A(N);


    Eigen::ConjugateGradient<ConstantSizeMatrix, Eigen::Lower | Eigen::Upper, Eigen::IdentityPreconditioner> cg;
    cg.compute(A);
    cg.setMaxIterations(maxIterations);
    cg.setTolerance(tolerance);

    Eigen::VectorXd u = b;
    {
        boost::timer::auto_cpu_timer t("cg took %w seconds\n");
        Eigen::ConjugateGradient<ConstantSizeMatrix, Eigen::Lower | Eigen::Upper, Eigen::IdentityPreconditioner> cg;
        cg.compute(A);
        u.noalias() = cg.solveWithGuess(b, u * 1.0 / 4.0);
        std::cout << "CG:       #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;
    }

    using MatrixT = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    std::cout << Eigen::Map<MatrixT>(u.data(), 2*N-1, N-1) << std::endl;



/*    pcl::PolygonMesh mesh = toMesh(u, N);
    std::cout << "i bims 1 lauchgemues" << std::endl;
    boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer (new pcl::visualization::PCLVisualizer ("3D Viewer"));

    viewer->setBackgroundColor (0, 0, 0);
    std::cout << "i bims wieder du honk" <<std::endl;
    viewer->addPolygonMesh(mesh,"meshes",0);
    std::cout << "i bims wieder 2" <<std::endl;
    viewer->addCoordinateSystem (1.0);
    viewer->initCameraParameters ();
    std::cout << "i bims wieder" <<std::endl;
    while (!viewer->wasStopped ())
    {
        viewer->spinOnce (100);
        boost::this_thread::sleep (boost::posix_time::microseconds (100000));
    }*/

}

