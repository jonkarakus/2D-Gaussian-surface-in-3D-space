using System;
using System.Collections.Generic;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Media;
using System.Windows.Media.Media3D;
using System.Windows.Media.Animation;

namespace Gaussian3D
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        private double sigma = 1.0;
        private double amplitude = 5.0;
        private int resolution = 150; // Higher resolution for better quality

        private Model3DGroup modelGroup;
        //private GeometryModel3D gaussianModel;
        private PerspectiveCamera camera;

        public MainWindow()
        {
            InitializeComponent();
            InitializeScene();
        }

        private void InitializeScene()
        {
            Viewport3D viewport = new Viewport3D();

            // Create camera and viewing angle
            camera = new PerspectiveCamera();
            camera.Position = new Point3D(15, 12, 15);
            camera.LookDirection = new Vector3D(-15, -12, -15);
            camera.UpDirection = new Vector3D(0, 1, 0);
            camera.FieldOfView = 45;
            viewport.Camera = camera;

            // Create model group to contain 3D objects
            modelGroup = new Model3DGroup();

            AddLighting();
            AddCoordinateSystem();
            AddGridPlane();
            CreateGaussianSurface();
            AddRotationAnimation();

            // Add model to the viewport
            ModelVisual3D modelVisual = new ModelVisual3D();
            modelVisual.Content = modelGroup;
            viewport.Children.Add(modelVisual);

            // Add viewport to the window
            MainViewportContainer.Children.Add(viewport);
        }

        private void AddLighting()
        {
            // Key light
            DirectionalLight keyLight = new DirectionalLight();
            keyLight.Color = Colors.White;
            keyLight.Direction = new Vector3D(-0.5, -0.5, -0.5);
            keyLight.Direction.Normalize();
            modelGroup.Children.Add(keyLight);

            // Fill light 
            DirectionalLight fillLight = new DirectionalLight();
            fillLight.Color = Color.FromRgb(220, 220, 255); // Slightly blue tint
            fillLight.Direction = new Vector3D(0.5, -0.2, 0.5);
            fillLight.Direction.Normalize();
            modelGroup.Children.Add(fillLight);

            // Rim light - highlights edges
            DirectionalLight rimLight = new DirectionalLight();
            rimLight.Color = Color.FromRgb(255, 255, 220); // Slightly yellow tint
            rimLight.Direction = new Vector3D(0.2, -0.8, -0.2);
            rimLight.Direction.Normalize();
            modelGroup.Children.Add(rimLight);

            // Ambient light 
            AmbientLight ambientLight = new AmbientLight();
            ambientLight.Color = Color.FromRgb(64, 64, 80); // Slightly blue ambient
            modelGroup.Children.Add(ambientLight);
        }

        private void AddCoordinateSystem()
        {
            double axisLength = 10;
            double axisThickness = 0.06;
            //XYZ Axis 
            modelGroup.Children.Add(CreateAxis(new Point3D(0, 0, 0), new Point3D(axisLength, 0, 0), Colors.Red, axisThickness));
            modelGroup.Children.Add(CreateAxis(new Point3D(0, 0, 0), new Point3D(0, axisLength, 0), Colors.Green, axisThickness));
            modelGroup.Children.Add(CreateAxis(new Point3D(0, 0, 0), new Point3D(0, 0, axisLength), Colors.Blue, axisThickness));
        }

        private void AddGridPlane()
        {
            // Create grid on the XZ plane (Y=0)
            double gridSize = 10;
            int gridLines = 20;
            double lineThickness = 0.02;

            for (int i = -gridLines / 2; i <= gridLines / 2; i++)
            {
                // X grid lines
                double pos = i * (gridSize / (gridLines / 2));
            
                    modelGroup.Children.Add(CreateAxis(
                        new Point3D(pos, 0, -gridSize),
                        new Point3D(pos, 0, gridSize),
                        Color.FromArgb(100, 180, 180, 180),
                        lineThickness));

                    modelGroup.Children.Add(CreateAxis(
                        new Point3D(-gridSize, 0, pos),
                        new Point3D(gridSize, 0, pos),
                        Color.FromArgb(100, 180, 180, 180),
                        lineThickness));
                
            }
        }

        private GeometryModel3D CreateAxis(Point3D start, Point3D end, Color color, double thickness)
        {
            Vector3D direction = end - start;
            double length = direction.Length;
            direction.Normalize();

            // Create cylinder along the specified axis
            MeshGeometry3D mesh = new MeshGeometry3D();
            int segments = 8;

            // Create perpendicular vector
            Vector3D perpendicular;
            if (Math.Abs(direction.X) < Math.Abs(direction.Y) && Math.Abs(direction.X) < Math.Abs(direction.Z))
                perpendicular = Vector3D.CrossProduct(new Vector3D(1, 0, 0), direction);
            else if (Math.Abs(direction.Y) < Math.Abs(direction.Z))
                perpendicular = Vector3D.CrossProduct(new Vector3D(0, 1, 0), direction);
            else
                perpendicular = Vector3D.CrossProduct(new Vector3D(0, 0, 1), direction);

            perpendicular.Normalize();
            Vector3D perpendicular2 = Vector3D.CrossProduct(direction, perpendicular);
            perpendicular2.Normalize();

            // circle of vertices at the start
            for (int i = 0; i < segments; i++)
            {
                double angle = 2 * Math.PI * i / segments;
                Vector3D offset = perpendicular * Math.Cos(angle) * thickness + perpendicular2 * Math.Sin(angle) * thickness;
                mesh.Positions.Add(start + offset);
            }

            // circle of vertices at the end
            for (int i = 0; i < segments; i++)
            {
                double angle = 2 * Math.PI * i / segments;
                Vector3D offset = perpendicular * Math.Cos(angle) * thickness + perpendicular2 * Math.Sin(angle) * thickness;
                mesh.Positions.Add(end + offset);
            }

            // Create triangles
            for (int i = 0; i < segments; i++)
            {
                int i0 = i;
                int i1 = (i + 1) % segments;
                int i2 = i + segments;
                int i3 = (i + 1) % segments + segments;

                mesh.TriangleIndices.Add(i0);
                mesh.TriangleIndices.Add(i1);
                mesh.TriangleIndices.Add(i2);

                mesh.TriangleIndices.Add(i2);
                mesh.TriangleIndices.Add(i1);
                mesh.TriangleIndices.Add(i3);
            }

            // Create material
            Material material = new DiffuseMaterial(new SolidColorBrush(color));

            // Create the model
            GeometryModel3D model = new GeometryModel3D(mesh, material);
            return model;
        }

        private void CreateGaussianSurface()
        {
            MeshGeometry3D mesh = new MeshGeometry3D();

            // Create gradient for material
            List<GeometryModel3D> surfaceModels = new List<GeometryModel3D>();

            // Parameters
            double halfSize = 6.0;
            double step = 2 * halfSize / resolution;

            // Create color map for Gaussian heights
            List<Color> colorMap = GenerateColorMap();

            // Create vertices for mesh
            Point3D[,] vertices = new Point3D[resolution + 1, resolution + 1];
            double[,] heights = new double[resolution + 1, resolution + 1];
            double maxHeight = 0;

            // Calculate all heights and find the maximum
            for (int i = 0; i <= resolution; i++)
            {
                double x = -halfSize + i * step;

                for (int j = 0; j <= resolution; j++)
                {
                    double z = -halfSize + j * step;
                    double y = GaussianFunction(x, z); // Note that Y is up in 3D space

                    vertices[i, j] = new Point3D(x, y, z);
                    heights[i, j] = y;
                    maxHeight = Math.Max(maxHeight, y);
                }
            }

            // Create triangles with height-based colors
            for (int i = 0; i < resolution; i++)
            {
                for (int j = 0; j < resolution; j++)
                {
                    // Get four corners of each quad
                    Point3D p00 = vertices[i, j];
                    Point3D p10 = vertices[i + 1, j];
                    Point3D p01 = vertices[i, j + 1];
                    Point3D p11 = vertices[i + 1, j + 1];

                    // Calculate average height for this quad
                    double avgHeight = (heights[i, j] + heights[i + 1, j] + heights[i, j + 1] + heights[i + 1, j + 1]) / 4.0;
                    double normalizedHeight = avgHeight / maxHeight;

                    Color quadColor = GetColorForHeight(normalizedHeight, colorMap);

                    MeshGeometry3D quadMesh = new MeshGeometry3D();

                    // Add four corners
                    quadMesh.Positions.Add(p00);
                    quadMesh.Positions.Add(p10);
                    quadMesh.Positions.Add(p01);
                    quadMesh.Positions.Add(p11);

                    // Add two triangles to form the quad
                    quadMesh.TriangleIndices.Add(0);
                    quadMesh.TriangleIndices.Add(1);
                    quadMesh.TriangleIndices.Add(2);

                    quadMesh.TriangleIndices.Add(2);
                    quadMesh.TriangleIndices.Add(1);
                    quadMesh.TriangleIndices.Add(3);

                    // Create normals for lighting
                    Vector3D normal1 = CalculateNormal(p00, p10, p01);
                    Vector3D normal2 = CalculateNormal(p01, p10, p11);

                    quadMesh.Normals.Add(normal1);
                    quadMesh.Normals.Add(normal1);
                    quadMesh.Normals.Add(normal1);
                    quadMesh.Normals.Add(normal2);

                    MaterialGroup materialGroup = new MaterialGroup();

                    DiffuseMaterial diffuseMat = new DiffuseMaterial(new SolidColorBrush(quadColor));
                    materialGroup.Children.Add(diffuseMat);

                    SpecularMaterial specMat = new SpecularMaterial(new SolidColorBrush(Colors.White), 50);
                    materialGroup.Children.Add(specMat);

                    // Create the model for this specific quad
                    GeometryModel3D quadModel = new GeometryModel3D(quadMesh, materialGroup);
                    surfaceModels.Add(quadModel);
                }
            }

            // Add all surface models to main model
            foreach (var model in surfaceModels)
            {
                modelGroup.Children.Add(model);
            }
        }

        private List<Color> GenerateColorMap()
        {
            // Color Gradient
            List<Color> colorMap = new List<Color>();

            colorMap.Add(Color.FromRgb(0, 0, 255));     // Blue
            colorMap.Add(Color.FromRgb(0, 128, 255));   // Light blue
            colorMap.Add(Color.FromRgb(0, 255, 255));   // Cyan
            colorMap.Add(Color.FromRgb(0, 255, 128));   // Teal
            colorMap.Add(Color.FromRgb(0, 255, 0));     // Green
            colorMap.Add(Color.FromRgb(128, 255, 0));   // Light green
            colorMap.Add(Color.FromRgb(255, 255, 0));   // Yellow
            colorMap.Add(Color.FromRgb(255, 128, 0));   // Orange
            colorMap.Add(Color.FromRgb(255, 0, 0));     // Red

            return colorMap;
        }

        private Color GetColorForHeight(double normalizedHeight, List<Color> colorMap)
        {
            if (normalizedHeight <= 0) return colorMap[0];
            if (normalizedHeight >= 1) return colorMap[colorMap.Count - 1];

            double scaledHeight = normalizedHeight * (colorMap.Count - 1);
            int lowerIndex = (int)Math.Floor(scaledHeight);
            int upperIndex = (int)Math.Ceiling(scaledHeight);

            if (lowerIndex == upperIndex) return colorMap[lowerIndex];

            double fraction = scaledHeight - lowerIndex;

            return InterpolateColors(colorMap[lowerIndex], colorMap[upperIndex], fraction);
        }

        private Color InterpolateColors(Color c1, Color c2, double fraction)
        {
            byte r = (byte)(c1.R + (c2.R - c1.R) * fraction);
            byte g = (byte)(c1.G + (c2.G - c1.G) * fraction);
            byte b = (byte)(c1.B + (c2.B - c1.B) * fraction);
            byte a = (byte)(c1.A + (c2.A - c1.A) * fraction);

            return Color.FromArgb(a, r, g, b);
        }

        private Vector3D CalculateNormal(Point3D p1, Point3D p2, Point3D p3)
        {
            Vector3D v1 = p2 - p1;
            Vector3D v2 = p3 - p1;

            Vector3D normal = Vector3D.CrossProduct(v1, v2);
            normal.Normalize();

            return normal;
        }

        private double GaussianFunction(double x, double z)
        {
            // Calculate standard Gaussian value
            double gaussianValue = amplitude * Math.Exp(-(x * x + z * z) / (2 * sigma * sigma));

            // Apply vertical offset
            double verticalOffset = amplitude * 0.2; 
            return gaussianValue - verticalOffset;
        }

        private void AddRotationAnimation()
        {
            System.Windows.Threading.DispatcherTimer timer = new System.Windows.Threading.DispatcherTimer();
            timer.Interval = TimeSpan.FromMilliseconds(16); // ~60 fps

            double angle = 0;
            double radius = Math.Sqrt(camera.Position.X * camera.Position.X + camera.Position.Z * camera.Position.Z);
            double height = camera.Position.Y;

            timer.Tick += (sender, e) => {
                // Update angle
                angle += 0.1; // Degrees per frame
                if (angle >= 360) angle = 0;

                // Calculate new camera position
                double radians = angle * Math.PI / 180;
                camera.Position = new Point3D(
                    radius * Math.Cos(radians),
                    height,
                    radius * Math.Sin(radians)
                );

                // Keep camera at center
                camera.LookDirection = new Vector3D(-camera.Position.X, -camera.Position.Y, -camera.Position.Z);
                camera.UpDirection = new Vector3D(0, 1, 0);
            };

            // Start timer
            timer.Start();
        }
    }
}
