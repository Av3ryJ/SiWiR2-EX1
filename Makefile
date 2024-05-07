CC = g++
CFLAGS = -O3 -Wall -Winline -Wshadow -std=c++17
LDFLAGS =

# Verzeichnis zum Speichern der ausführbaren Dateien
BINDIR = bin

# Quelldateien
SOURCES = mgsolve.cpp grid.cpp
# Ersetzen der Dateiendungen .cpp durch .o für Objektdateien
OBJECTS = $(SOURCES:.cpp=.o)
# Name des ausführbaren Programms
EXECUTABLE = $(BINDIR)/mgsolve

# Standard-Ziel, das durch Ausführen von "make" erstellt wird
all: $(EXECUTABLE)

# Regel zum Erstellen des ausführbaren Programms
$(EXECUTABLE): $(OBJECTS)
	@mkdir -p $(BINDIR)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

# Regel zum Kompilieren der Quelldateien
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

# Aufräumen der Objektdateien und des ausführbaren Programms
clean:
	rm -f $(OBJECTS) $(EXECUTABLE)