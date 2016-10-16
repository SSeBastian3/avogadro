/******************************************************************************

  This source file is part of the Avogadro project.

  Copyright 2016 Kitware, Inc.

  This source code is released under the New BSD License, (the "License").

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

******************************************************************************/

#ifndef YAEHMOP_EXTENSION_H
#define YAEHMOP_EXTENSION_H

#include <avogadro/extension.h>

#include <QMutex>

namespace Avogadro {

  class YaehmopExtension : public Extension
  {
    Q_OBJECT
    AVOGADRO_EXTENSION("Yaehmop", tr("Yaehmop"),
                       tr("Yaehmop extension"))

  public:
    YaehmopExtension(QObject *parent=0);
    virtual ~YaehmopExtension();

    virtual QList<QAction *> actions() const;
    virtual QString menuPath(QAction *action) const;

    virtual void setMolecule(Molecule *molecule);

    virtual QDockWidget * dockWidget();
    virtual QUndoCommand* performAction(QAction *action, GLWidget *widget);

    // Write to Yaehmop and receive the output. This call is blocking (it will
    // not return until Yaehmop finishes). Returns true on success and false on
    // failure. The output is always written to QString output.
    bool executeYaehmop(QString input, QString& output) const;

  public slots:
    void calculateBandStructure() const;
    void calculateTotalDOS() const;
    void plotPartialDOS() const;
    void setParametersFile();
    void executeCustomInput() const;

  private:
    // @param displayBandData This will be set to true if we are to
    //                        display the band data for the user.
    // @param limitY Should we limit the y range?
    // @param minY MinY if we are limiting the y range.
    // @param maxY MaxY if we are limiting the y range.
    // @return The band calculation input.
    QString createYaehmopBandInput(bool& displayBandData, bool& limitY,
                                   double& minY, double& maxY) const;

    // @param displayDOSData This will be set to true if we are to
    //                       display the DOS data for the user.
    // @param useSmoothing This will be set to true if we are to
    //                     use Gaussian smoothing on the data.
    // @param stepE If useSmoothing is true, this will contain
    //              the energy step size to be used for Gaussian smoothing.
    // @param broadening If useSmoothing is true, this will contain
    //                   the broadening to be used for Gaussian smoothing.
    // @param limitY Should we limit the y-range?
    // @param minY MinY if we are limiting the y-range.
    // @param maxY MaxY if we are limiting the y-range.
    // @return The total DOS calculation input.
    QString createYaehmopTotalDOSInput(bool& displayDOSData,
                                       bool& useSmoothing,
                                       double& stepE,
                                       double& broadening, bool& limitY,
                                       double& minY, double& maxY) const;

    QString createGeometryAndLatticeInput() const;

    // Assuming a constant x difference, integrate the data using the
    // trapezoid rule and return it as a QList.
    // Each point adds the last point to it. The final value
    // in the QList is the total integration.
    // @param xDiff The distance between x values.
    // @param y The y value data.
    // @return The integration data.
    static QList<double> integrateDataTrapezoidal(double xDist,
                                                  const QVector<double>& y);

    // Smooths data using Gaussian smoothing.
    // @param densities The x values to be smoothed
    // @param energies The y values to be smoothed
    // @param stepE The new distance between energy points
    // @param broadening The broadening for the smoothing
    static void smoothData(QVector<double>& densities,
                           QVector<double>& energies,
                           double stepE, double broadening);

 public:
    // All of the following functions are just for saving settings during
    // while the program is running. All settings get reset when the program
    // is closed.
    // Perhaps a mutex is not necessary, but I added it for safety anyways...
    static void setBandNumKPoints(size_t numKPoints)
    {
      lock();
      s_bandNumKPoints = numKPoints;
      unlock();
    };

    static size_t getBandNumKPoints()
    {
      lock();
      size_t ret = s_bandNumKPoints;
      unlock();
      return ret;
    }

    static void setDOSKPoints(QString kpoints)
    {
      lock();
      s_dosKPoints = kpoints;
      unlock();
    };

    static QString getDOSKPoints()
    {
      lock();
      QString ret = s_dosKPoints;
      unlock();
      return ret;
    }

    static void setUseBroadening(bool useBroadening)
    {
      lock();
      s_useBroadening = useBroadening;
      unlock();
    };

    static bool getUseBroadening()
    {
      lock();
      bool ret = s_useBroadening;
      unlock();
      return ret;
    }

    static void setEnergyStepSize(double eStep)
    {
      lock();
      s_energyStepSize = eStep;
      unlock();
    };

    static double getEnergyStepSize()
    {
      lock();
      double ret = s_energyStepSize;
      unlock();
      return ret;
    }

    static void setBroadening(double broadening)
    {
      lock();
      s_broadening = broadening;
      unlock();
    };

    static double getBroadening()
    {
      lock();
      double ret = s_broadening;
      unlock();
      return ret;
    }

    static void setDisplayData(bool displayData)
    {
      lock();
      s_displayData = displayData;
      unlock();
    };

    static bool getDisplayData()
    {
      lock();
      bool ret = s_displayData;
      unlock();
      return ret;
    }

    static void setLimitY(bool limitY)
    {
      lock();
      s_limitY = limitY;
      unlock();
    };

    static bool getLimitY()
    {
      lock();
      bool ret = s_limitY;
      unlock();
      return ret;
    }

    static void setMinY(double minY)
    {
      lock();
      s_minY = minY;
      unlock();
    };

    static double getMinY()
    {
      lock();
      double ret = s_minY;
      unlock();
      return ret;
    }

    static void setMaxY(double maxY)
    {
      lock();
      s_maxY = maxY;
      unlock();
    };

    static double getMaxY()
    {
      lock();
      double ret = s_maxY;
      unlock();
      return ret;
    }

 private:

    static void lock() { s_mutex.lock(); };
    static void unlock() { s_mutex.unlock(); };

    static QMutex s_mutex;

    static size_t s_bandNumKPoints;
    static QString s_dosKPoints;
    static bool s_useBroadening;
    static double s_energyStepSize;
    static double s_broadening;
    static bool s_displayData;
    static bool s_limitY;
    static double s_minY;
    static double s_maxY;

    QList<QAction *> m_actions;
    Molecule *m_molecule;
    QString m_parametersFile;

  private Q_SLOTS:

  };

  class YaehmopExtensionFactory : public QObject, public PluginFactory
  {
    Q_OBJECT
    Q_INTERFACES(Avogadro::PluginFactory)
        AVOGADRO_EXTENSION_FACTORY(YaehmopExtension)
  };

} // end namespace Avogadro

#endif // YAEHMOP_EXTENSION_H
