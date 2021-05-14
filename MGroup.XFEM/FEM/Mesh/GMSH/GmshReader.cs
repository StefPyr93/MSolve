using System;
using System.Collections.Generic;
using System.IO;

//TODO: Enforce contiguous vertex ids
namespace MGroup.XFEM.FEM.Mesh.GMSH
{
    /// <summary>
    /// Creates meshes by reading GMSH output files (.msh). Unrecognized GMSH cell types will be ignored when parsing the 
    /// .msh file, therefore some care is needed.
    /// </summary>
    public class GmshReader: IMeshGenerator, IDisposable
    {
        private readonly StreamReader reader;
        private readonly int minDimensionOfCells;

        /// <summary>
        /// Opens the .msh file but doesn't read it.
        /// </summary>
        /// <param name="mshFilePath">
        /// The absolute path of the .msh file where GMSH has written the mesh data. The .msh file will not be modified.
        /// </param>
        /// <param name="minDimensionOfCells">
        /// The minimum dimension of cells that will be kept. E.g. if a 3D mesh is read and 
        /// <paramref name="minDimensionOfCells"/>=2, then 3D cells (bulk elements) and 2D cells (face cells of geometric 
        /// entities) will be kept, while 1D cells (edge cells of geometric entities) will be ignored.
        /// </param>
        public GmshReader(string mshFilePath, int minDimensionOfCells)
        {
            reader = new StreamReader(mshFilePath);

            if ((minDimensionOfCells != 1) && (minDimensionOfCells != 2) && (minDimensionOfCells != 3))
            {
                throw new ArgumentException(
                    "The min dimension of cells that will be kept must be 1, 2 or 3, but was " + minDimensionOfCells);
            }
            this.minDimensionOfCells = minDimensionOfCells;
        }

        ~GmshReader()
        {
            if (reader != null) reader.Dispose();
        }

        /// <summary>
        /// Reads the whole .msh file and converts it to MSolve mesh data.
        /// </summary>
        public PreprocessingMesh CreateMesh()
        {
            // Vertices must be listed before cells
            var mesh = new PreprocessingMesh();
            ReadVertices(mesh);
            ReadCells(mesh);
            return mesh;
        }


        public void Dispose()
        {
            if (reader != null) reader.Dispose();
        }

        private void ReadVertices(PreprocessingMesh mesh)
        {
            string line;

            // Find node segment
            while (true)
            {
                line = reader.ReadLine();
                if (line[0] == '$')
                    if (line.Equals("$Nodes")) break; // Next line will be the nodes count.
            }

            // Read the vertices
            int numVertices = int.Parse(reader.ReadLine());
            mesh.Vertices = new List<Vertex>(numVertices);
            while (true)
            {
                line = reader.ReadLine();
                if (line[0] == '$') break; // This line is "$EndNodes". Next line will be the next segment.
                else
                {
                    // Line: nodeID x y z
                    string[] words = line.Split(new char[] { ' ' });
                    int id = int.Parse(words[0]) - 1; // MSolve uses 0-based indexing
                    if (id != mesh.Vertices.Count)
                    {
                        throw new NotImplementedException("GMSH vertex ids must be contiguous");
                    }
                    mesh.Vertices.Add(new Vertex(id, double.Parse(words[1]), double.Parse(words[2]), double.Parse(words[3])));
                }
            }
        }

        // It must be called after vertices are read.
        private void ReadCells(PreprocessingMesh mesh)
        {
            string line;

            // Find cell segment of the file
            while (true)
            {
                line = reader.ReadLine();
                if (line[0] == '$')
                    if (line.Equals("$Elements")) break; // Next line will be the number of cells.
            }

            // Decide which cell types will be kept
            var cellFactory = new GmshCellFactory(mesh.Vertices);
            HashSet<int> validCellTypes = cellFactory.GetValidGmshCellTypes(minDimensionOfCells);

            // Read the cells
            int numTotalCells = int.Parse(reader.ReadLine()); // not all of them will be kept
            mesh.Cells = new List<Cell>(numTotalCells);
            while (true)
            {
                line = reader.ReadLine();
                if (line[0] == '$') break; // This line is "$EndElements". Next line will be the next segment.
                else
                {
                    // Line: cellID cellType tagsCount <tags>(0 or more) vertexIds(2 or more)
                    string[] words = line.Split(new char[] { ' ' });
                    int gmshCellType = int.Parse(words[1]);
                    if (!validCellTypes.Contains(gmshCellType)) continue; // Ignore unrecognized cells, instead of throwing an exception

                    int numTags = int.Parse(words[2]);
                    int physicalGroup = int.Parse(words[3]);
                    int firstVertexPos = 3 + numTags;
                    int numVertices = words.Length - firstVertexPos;

                    int[] vertexIDs = new int[numVertices];
                    for (int i = 0; i < numVertices; ++i)
                    {
                        vertexIDs[i] = int.Parse(words[firstVertexPos + i]) - 1; // MSolve uses 0-based indexing
                    }

                    // Use contiguous numbering for cell IDs, instead of the one provided by GMSH.
                    Cell cell = cellFactory.CreateGmshCell(mesh.Cells.Count, gmshCellType, vertexIDs, physicalGroup);
                    mesh.Cells.Add(cell);
                }
            }
        }
    }
}
