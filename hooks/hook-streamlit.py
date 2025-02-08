from PyInstaller.utils.hooks import copy_metadata

datas = []

# Include metadata for all required packages
datas += copy_metadata("streamlit")
datas += copy_metadata("extras")
datas += copy_metadata("numpy")
datas += copy_metadata("pandas_flavor")
datas += copy_metadata("plotly")
datas += copy_metadata("openpyxl")
datas += copy_metadata("scipy") 
datas += copy_metadata("networkx")
datas += copy_metadata("statsmodels")



