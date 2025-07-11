python3 -m venv venv
source venv/bin/activate
pip install pyvista scipy numpy
/Applications/MATLAB_R2024a.app/bin/matlab -nodisplay -nosplash -r "run_turbulence; exit"
deactivate
rm -rf venv
