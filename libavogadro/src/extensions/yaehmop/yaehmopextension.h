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
    // @return The band calculation input.
    QString createYaehmopBandInput(bool& displayBandData) const;

    // @param displayDOSData This will be set to true if we are to
    //                       display the DOS data for the user.
    // @param useSmoothing This will be set to true if we are to
    //                     use Gaussian smoothing on the data.
    // @param stepE If useSmoothing is true, this will contain
    //              the energy step size to be used for Gaussian smoothing.
    // @param broadening If useSmoothing is true, this will contain
    //                   the broadening to be used for Gaussian smoothing.
    // @return The total DOS calculation input.
    QString createYaehmopTotalDOSInput(bool& displayDOSData,
                                       bool& useSmoothing,
                                       double& stepE,
                                       double& broadening) const;

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
