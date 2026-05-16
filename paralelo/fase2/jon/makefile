# 1. Variables de compilador y banderas
CC       = mpicc
CFLAGS   = -Wall -O2
TARGET   = ejecutable

# 2. Lista de archivos fuente y objetos
SRCS     = faux_p.c heat_fase2_p.c diffusion_fase2_p.c
OBJS     = $(SRCS:.c=.o)

# 3. Regla principal (la que se ejecuta con solo escribir 'make')
all: $(TARGET)

# 4. Cómo enlazar los objetos para crear el ejecutable
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) -lm

# 5. Cómo compilar cada archivo .c en un .o
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# 6. eliminar archivos .o
    clean_objs:
	rm -f $(OBJS)

# 7. Limpieza de archivos temporales
clean:
	rm -f $(OBJS) $(TARGET)
