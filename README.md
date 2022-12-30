# met_antartica
Scripts elaborados para o trabalho final da disciplina meteorologia antártica

# Observações: 

Os scripts disponibilizados foram escritos para serem utilizados em uma máquina com sistema operacional windows. Uma vez que há diferenças entre os sistemas, segue um exemplo de como modificar o código para que não hajam problemas ao utiliza-los.

# 1) windows (forma escrita nos scripts)

file_1 = xr.open_dataset(r'C:\Users\User\OneDrive\Documentos\Meteorologia\bia\antartica\dados\antartica-pressure-levels.nc').metpy.parse_cf()
    
# 2) linux/ubuntu

file_1 = xr.open_dataset('/Users/User/OneDrive/Documentos/Meteorologia/bia/antartica/dados/antartica-pressure-levels.nc').metpy.parse_cf()
