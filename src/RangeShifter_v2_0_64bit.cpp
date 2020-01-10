//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop
#include <tchar.h>
//---------------------------------------------------------------------------
USEFORM("VCL\FormMain.cpp", frmMain);
USEFORM("VCL\FormLand.cpp", frmLand);
USEFORM("VCL\FormSeeding.cpp", frmSeeding);
USEFORM("VCL\FormMove.cpp", frmMove);
USEFORM("VCL\FormGenetics.cpp", frmGenetics);
USEFORM("VCL\FormDensity.cpp", frmDensity);
USEFORM("VCL\FormArtificialLand.cpp", frmGenerateLand);
USEFORM("VCL\FormEnvGradient.cpp", frmEnvGradient);
USEFORM("VCL\FormDynLand.cpp", frmDynLand);
USEFORM("VCL\FormVisualTraits1.cpp", frmVisualTraits1);
USEFORM("VCL\FormVisualTraits0.cpp", frmVisualTraits0);
USEFORM("VCL\FormVisualTraits2.cpp", frmVisualTraits2);
USEFORM("VCL\FormVisualPatch.cpp", frmVisualPatch);
USEFORM("VCL\FormSpecies.cpp", frmSpecies);
USEFORM("VCL\FormSim.cpp", frmSim);
USEFORM("VCL\FormVisualGrad.cpp", frmVisualGrad);
USEFORM("VCL\FormVisualCost.cpp", frmVisualCost);
//---------------------------------------------------------------------------
int WINAPI _tWinMain(HINSTANCE, HINSTANCE, LPTSTR, int)
{
	try
	{
		Application->Initialize();
		Application->MainFormOnTaskBar = true;
		Application->Title = "RangeShifter_v2.0";
		Application->CreateForm(__classid(TfrmMain), &frmMain);
		Application->CreateForm(__classid(TfrmGenerateLand), &frmGenerateLand);
		Application->CreateForm(__classid(TfrmDensity), &frmDensity);
		Application->CreateForm(__classid(TfrmEnvGradient), &frmEnvGradient);
		Application->CreateForm(__classid(TfrmLand), &frmLand);
		Application->CreateForm(__classid(TfrmMove), &frmMove);
		Application->CreateForm(__classid(TfrmSeeding), &frmSeeding);
		Application->CreateForm(__classid(TfrmSim), &frmSim);
		Application->CreateForm(__classid(TfrmSpecies), &frmSpecies);
		Application->CreateForm(__classid(TfrmVisualCost), &frmVisualCost);
		Application->CreateForm(__classid(TfrmVisualGrad), &frmVisualGrad);
		Application->CreateForm(__classid(TfrmVisualPatch), &frmVisualPatch);
		Application->CreateForm(__classid(TfrmDynLand), &frmDynLand);
		Application->CreateForm(__classid(TfrmGenetics), &frmGenetics);
		Application->CreateForm(__classid(TfrmVisualTraits0), &frmVisualTraits0);
		Application->CreateForm(__classid(TfrmVisualTraits1), &frmVisualTraits1);
		Application->CreateForm(__classid(TfrmVisualTraits2), &frmVisualTraits2);
		Application->Run();
	}
	catch (Exception &exception)
	{
		Application->ShowException(&exception);
	}
	catch (...)
	{
		try
		{
			throw Exception("");
		}
		catch (Exception &exception)
		{
			Application->ShowException(&exception);
		}
	}
	return 0;
}
//---------------------------------------------------------------------------
