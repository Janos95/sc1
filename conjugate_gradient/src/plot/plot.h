#include <vtkVersion.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>
#include <vtkSmartPointer.h>
#include <vtkChartXY.h>
#include <vtkTable.h>
#include <vtkPlot.h>
#include <vtkFloatArray.h>
#include <vtkContextView.h>
#include <vtkContextScene.h>
#include <vtkPen.h>

/*constexpr red = 
constexpr b*/l


class Plotter
{
public:
    
    Plotter()
    {
        view = vtkSmartPointer<vtkContextView>::New();
        view->GetRenderer()->SetBackground(1.0, 1.0, 1.0);

        chart = vtkSmartPointer<vtkChartXY>::New();
        view->GetScene()->AddItem(chart);
    }

    template<typename Range>
    void setXAxis(const Range range)
    {
        for(const auto& x: range)
            
    }
    
    template<typename F>
    void plot(F f, const std::array<char>& color, )
    {

        vtkPlot *line = chart->AddPlot(vtkChart::LINE);

        vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();



        line->SetInputData(table, 0, 1);
        line->SetColor(color[0], color[2], color[3], 255);
        line->SetWidth(1.0);

    }
        
    void spin()
    { 
        view->GetInteractor()->Initialize();
        view->GetInteractor()->Start();
    }
    
private:
    vtkSmartPointer<vtkContextView> view;
    vtkSmartPointer<vtkChartXY> chart
    
    std::vector<float> xAxis;
    std::vector<std::vector<float>> ys;
};









int main(int, char *[])
{
  vtkSmartPointer<vtkTable> table = 
    vtkSmartPointer<vtkTable>::New();

  vtkSmartPointer<vtkFloatArray> arrX = 
    vtkSmartPointer<vtkFloatArray>::New();
  arrX->SetName("X Axis");
  table->AddColumn(arrX);

  vtkSmartPointer<vtkFloatArray> arrC = 
    vtkSmartPointer<vtkFloatArray>::New();
  arrC->SetName("Cosine");
  table->AddColumn(arrC);

  vtkSmartPointer<vtkFloatArray> arrS = 
    vtkSmartPointer<vtkFloatArray>::New();
  arrS->SetName("Sine");
  table->AddColumn(arrS);

  // Fill in the table with some example values
  int numPoints = 69;
  float inc = 7.5 / (numPoints-1);
  table->SetNumberOfRows(numPoints);
  for (int i = 0; i < numPoints; ++i)
  {
    table->SetValue(i, 0, i * inc);
    table->SetValue(i, 1, cos(i * inc));
    table->SetValue(i, 2, sin(i * inc));
  }

  // Set up the view
//   vtkSmartPointer<vtkContextView> view = 
//     vtkSmartPointer<vtkContextView>::New();
//   view->GetRenderer()->SetBackground(1.0, 1.0, 1.0);

  // Add multiple line plots, setting the colors etc
 // Add multiple line plots, setting the colors etc

  // (ifdef-ed out on Windows because DASH_LINE does not work on Windows
  //  machines with built-in Intel HD graphics card...)

  //view->GetRenderWindow()->SetMultiSamples(0);

  // Start interactor
  view->GetInteractor()->Initialize();
  view->GetInteractor()->Start();

  return EXIT_SUCCESS;
