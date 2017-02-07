TEMPLATE = app
TARGET = 
DEPENDPATH += .
INCLUDEPATH += "/usr/local/include"
LIBS += -L"/usr/lib" -lOpenCL \
	

QT           += opengl
CONFIG       += console
#RESOURCES     = Ipm.qrc
HEADERS += io/*.h io/LibJpeg/*.h imageModel/Point.h imageModel/Line.h imageModel/Edge.h imageModel/Matrix.h imageModel/Image.h \
			segmentation/Thresholds.h segmentation/Canny.h segmentation/Suzuki.h histograms/*.h pht/*.h \
		  pointInterest/*.h correlation/*.h utils/*.h ui/*.h MAELab.h 
	
SOURCES += MAELab_CI.cpp io/*.cpp io/LibJpeg/*.c imageModel/Point.cpp imageModel/Line.cpp imageModel/Edge.cpp imageModel/Matrix.cpp imageModel/Image.cpp \
		segmentation/Thresholds.cpp segmentation/Canny.cpp segmentation/Suzuki.cpp \
		 histograms/*.cpp pht/*.cpp pointInterest/*.cpp correlation/*.cpp \
		 utils/*.cpp ui/*.cpp MAELab.cpp 