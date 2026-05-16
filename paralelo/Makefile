# --- Variables originales tuyas ---
CC       = mpicc
CFLAGS   = -Wall -O2
TARGET   = bin

# --- Phony targets ---
.PHONY: clean FORCE

# --- REGLA SALVAVIDAS ---
# Le decimos a Make explícitamente que no intente compilar el propio Makefile
Makefile: ;

# --- Regla de limpieza ---
clean:
	@echo "Limpiando ejecutables y restos..."
	@rm -f *.o
	@rm -rf $(TARGET)
	@echo "¡Directorio limpio!"

# --- Regla mágica (Comodín) ---
%: FORCE
	@# 1. Limpiamos la entrada: quitamos el "./" del principio y la "/" del final
	$(eval RAW_DIR = $(patsubst ./%,%,$@))
	$(eval DIR = $(patsubst %/,%,$(RAW_DIR)))

	@# 2. Asignamos TUS archivos fuente usando el directorio limpio
	$(eval SRCS = $(DIR)/faux.c $(DIR)/heat.c $(DIR)/diffusion.c)

	@# 3. Generamos el nombre del ejecutable (ej: fase2 o fase1-asincrono)
	$(eval EXEC = $(subst /,-,$(DIR)))

	@# 4. Comprobamos que el directorio y el código existan
	@if [ ! -f "$(DIR)/faux.c" ]; then \
		echo "Error: No se encuentra $(DIR)/faux.c. ¿Seguro que es el directorio correcto?"; \
	else \
		echo "--> 1. Compilando archivos objeto de $(DIR) con mpicc..."; \
		$(CC) $(CFLAGS) -c $(SRCS); \
		echo "--> 2. Creando carpeta '$(TARGET)' y enlazando el ejecutable '$(EXEC)'..."; \
		mkdir -p $(TARGET); \
		$(CC) *.o -o $(TARGET)/$(EXEC); \
		echo "--> 3. Eliminando archivos objeto (.o) automáticamente..."; \
		rm -f *.o; \
		echo "==> ¡COMPILACIÓN EXITOSA! Ejecutable listo en: $(TARGET)/$(EXEC)"; \
	fi

# Obliga a Make a evaluar la regla siempre
FORCE:
