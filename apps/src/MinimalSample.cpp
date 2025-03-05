/*

Copyright (c) 2005-2023, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

/**
 * @file
 *
 * This file gives an example of how you can create your own executable
 * in a user project.
 */

#include <iostream>
#include <string>

// Must be included before any other serialization headers
#include "CheckpointArchiveTypes.hpp"

#include "ArchiveLocationInfo.hpp"
#include "ArchiveOpener.hpp"
#include "FileFinder.hpp"
#include "SimulationTime.hpp"

#include "Exception.hpp"
#include "ExecutableSupport.hpp"
#include "PetscException.hpp"
#include "PetscTools.hpp"

#include "CellBasedSimulationArchiver.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "DiffusionForce.hpp"
#include "ElongationTrackingModifier.hpp"
#include "FarhadifarForce.hpp"
#include "NoCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"
#include "RandomCellKiller.hpp"
#include "RandomNumberGenerator.hpp"
#include "SmartPointers.hpp"
#include "StemCellProliferativeType.hpp"
#include "Toroidal2dVertexMesh.hpp"
#include "TransitCellProliferativeType.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "VoronoiVertexMeshGenerator.hpp"


int main(int argc, char* argv[])
{
    // This sets up PETSc and prints out copyright information, etc.
    // ExecutableSupport::StandardStartup(&argc, &argv);
    ExecutableSupport::StartupWithoutShowingCopyright(&argc, &argv);

    int exit_code = ExecutableSupport::EXIT_OK;

    // You should put all the main code within a try-catch, to ensure that
    // you clean up PETSc before quitting.
    try
    {
        double A0 = 0.85; // target area

        // Parameters for initialization
        size_t M = 3; // Number of cells in confluent tissue to be simulated

        long long random_seed = -1;

        // Simulation step parameters
        double total_production_run_time = 0.5;
        double simulation_delta_step = 0.002;
        double simulation_sampling_timestep = 0.5;

        // Input/Output parameters
        std::string output_directory = "./output_vm/"; // Path prefix for output files to be produced

        RandomNumberGenerator* random_generator = RandomNumberGenerator::Instance();
        random_generator->Reseed(random_seed);

        size_t simulation_iterations_per_sampling = roundf(simulation_sampling_timestep / simulation_delta_step);

        std::cerr << "Setting output frequency to " << simulation_iterations_per_sampling << " steps" << std::endl;

        // TODO: Deal with confluent and finite simulations differently
        if (PetscTools::AmMaster())
        {
            SimulationTime::Instance()->SetStartTime(0.0);
            VertexBasedCellPopulation<2>* cell_population = nullptr;
            OffLatticeSimulation<2>* simulator = nullptr;
            boost::shared_ptr<Toroidal2dVertexMesh> p_mesh;

            MAKE_PTR(ElongationTrackingModifier<2>, p_modifier_elongation);

            VoronoiVertexMeshGenerator generator(M, M, 0, A0); // Parameters are: cells across, cells up
            p_mesh = generator.GetToroidalMesh();

            MAKE_PTR(StemCellProliferativeType, p_stem_type);
            std::vector<CellPtr> cells;
            CellsGenerator<NoCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_stem_type);

            cell_population = new VertexBasedCellPopulation<2>(*p_mesh, cells);

            // Calculate ESF
            p_modifier_elongation->UpdateCellData(*cell_population);

            simulator = new OffLatticeSimulation<2>(*cell_population);
            simulator->SetOutputDirectory(output_directory);

            simulator->SetSamplingTimestepMultiple(simulation_iterations_per_sampling);
            simulator->SetDt(simulation_delta_step);

            CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(simulator);


            if (total_production_run_time > 0.0)
            {
                simulator->SetEndTime(total_production_run_time);
                simulator->Solve();
            }

            // Calculate ESF
            p_modifier_elongation->UpdateCellData(*cell_population);

            std::cerr << "Saving state..." << std::endl;
            // Classical saving
            CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(simulator);
            std::cerr << "Finished saving of final state..." << std::endl;

            delete simulator;
            simulator = nullptr;

            // Classical loading
            simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load(output_directory, total_production_run_time);
            std::cerr << "Finished classical loading" << std::endl;

            std::cerr << "Cleaning up ..." << std::endl;

            SimulationTime::Instance()->Destroy();
            delete simulator;
            if (cell_population != nullptr)
            {
                delete cell_population;
            }
        }
    }
    catch (const Exception& e)
    {
        ExecutableSupport::PrintError(e.GetMessage());
        exit_code = ExecutableSupport::EXIT_ERROR;
    }

    // Optional - write the machine info to file.
    ExecutableSupport::WriteMachineInfoFile("machine_info");

    // End by finalizing PETSc, and returning a suitable exit code.
    // 0 means 'no error'
    ExecutableSupport::FinalizePetsc();
    return exit_code;
}
