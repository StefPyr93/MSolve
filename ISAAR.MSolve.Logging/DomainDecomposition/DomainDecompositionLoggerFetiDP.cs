using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.Augmentation;

namespace ISAAR.MSolve.Logging.DomainDecomposition
{
    public class DomainDecompositionLoggerFetiDP : IDomainDecompositionLogger
    {
        private readonly string plotDirectoryPath;
        private readonly ICornerNodeSelection cornerNodesSelection;
        private readonly IMidsideNodesSelection midsideNodesSelection;
        private readonly bool shuffleSubdomainColors;
        private int analysisStep;

        //TODO: make sure the path does not end in "\"
        public DomainDecompositionLoggerFetiDP(string plotDirectoryPath, ICornerNodeSelection cornerNodesSelection, 
            IMidsideNodesSelection midsideNodesSelection = null, bool shuffleSubdomainColors = false) 
        {
            this.plotDirectoryPath = plotDirectoryPath;
            this.cornerNodesSelection = cornerNodesSelection;
            this.midsideNodesSelection = midsideNodesSelection;
            this.shuffleSubdomainColors = shuffleSubdomainColors;
            analysisStep = 0;
        }

        public void PlotSubdomains(IModel model)
        {
            var writer = new MeshPartitionWriter(shuffleSubdomainColors);
            writer.WriteSubdomainElements($"{plotDirectoryPath}\\subdomains_{analysisStep}.vtk", model);
            writer.WriteBoundaryNodes($"{plotDirectoryPath}\\boundary_nodes_{analysisStep}.vtk", model);

            INode[] crosspoints = model.EnumerateNodes().Where(n => n.SubdomainsDictionary.Count > 2).ToArray();
            writer.WriteSpecialNodes($"{plotDirectoryPath}\\crosspoints_{analysisStep}.vtk", "crosspoints", crosspoints);

            writer.WriteSpecialNodes($"{plotDirectoryPath}\\corner_nodes_{analysisStep}.vtk", "corner_nodes", 
                cornerNodesSelection.GlobalCornerNodes);

            if (midsideNodesSelection != null)
            {
                writer.WriteSpecialNodes($"{plotDirectoryPath}\\midside_nodes_{analysisStep}.vtk", "midside_nodes",
                midsideNodesSelection.MidsideNodesGlobal);
            }

            INode[] constrainedNodes = model.EnumerateNodes().Where(node => node.Constraints.Count > 0).ToArray();
            writer.WriteSpecialNodes($"{plotDirectoryPath}\\constrained_nodes_{analysisStep}.vtk", "constrained_nodes",
                constrainedNodes);

            ++analysisStep;
        }
    }
}
