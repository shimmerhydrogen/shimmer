#include <iostream>
#include <QMenuBar>
#include <QStatusBar>
#include <QFileDialog>

#include "CShimmerMainWindow.h"
#include "infra/infrastructure.h"
#include "errors.h"

CShimmerCentralWidget::CShimmerCentralWidget(QWidget *parent = nullptr)
    : QWidget(parent)
{

}


CShimmerMainWindow::CShimmerMainWindow(QWidget *parent)
{
    fileMenu_ = new QMenu("&File");
    QAction* quitaction = fileMenu_->addAction(tr("&Quit"));
    menuBar()->addMenu(fileMenu_);

    simMenu_ = new QMenu("&Simulation");
    QAction* simrunaction = simMenu_->addAction(tr("&Load and run"));
    menuBar()->addMenu(simMenu_);

    dbMenu_ = new QMenu("&Database");
    QAction* initdbaction = dbMenu_->addAction(tr("&Initialize new"));
    QAction* makedemodbaction = dbMenu_->addAction(tr("&Make demo DB"));
    QAction* testdbaction = dbMenu_->addAction(tr("&Load and test"));
    menuBar()->addMenu(dbMenu_);

    setCentralWidget( new CShimmerCentralWidget() );

    QObject::connect(quitaction, SIGNAL(triggered()), this, SLOT(qa()));
    QObject::connect(initdbaction, SIGNAL(triggered()), this, SLOT(database_init()));
    QObject::connect(testdbaction, SIGNAL(triggered()), this, SLOT(database_test()));
    QObject::connect(simrunaction, SIGNAL(triggered()), this, SLOT(simulation_run()));
    statusBar()->showMessage(tr("Ready"));
}

void
CShimmerMainWindow::database_init(void)
{
    const size_t BUFSIZE = 1024;
    QFileDialog fileDialog;
    fileDialog.setDefaultSuffix(QString("db"));
    fileDialog.setAcceptMode(QFileDialog::AcceptSave);
    fileDialog.setFileMode(QFileDialog::AnyFile);
    fileDialog.setNameFilter(tr("SQLite3 DB (*.db)"));
    

    if (not fileDialog.exec())
        return;

    QStringList filenames = fileDialog.selectedFiles();
    if (filenames.size() == 0)
        return;
    
    std::string filename =
        filenames[ filenames.size()-1 ].toStdString();

    char buf[BUFSIZE];
    FILE *fh = nullptr;
    FILE *ph = nullptr;

    if (filename.length() == 0) {
        statusBar()->showMessage("Invalid file name");
        return;
    }

    if ( (fh = fopen("../../sqlite/shimmer.sql", "r")) )
        goto foundok;

    if ( (fh = fopen("shimmer.sql", "r")) )
        goto foundok;

        statusBar()->showMessage("Failed to open SQL schema. Corrupt installation?");
    return;
    
foundok:
    unlink(filename.c_str());
    std::string popenstr = "sqlite3 " + filename;

    ph = popen(popenstr.c_str(), "w");
    if (not ph) {
        statusBar()->showMessage("Failed to populate DB: system error");
        return;
    }

    while (!feof(fh)) {
        size_t nread = fread(buf, 1, BUFSIZE, fh);
        fwrite(buf, 1, nread, ph);
    }

    pclose(fh);
    int sqlite_ret = pclose(ph);
    if (sqlite_ret != 0) {
        statusBar()->showMessage("Failed to populate DB: SQLite error");
        return;
    }


    statusBar()->showMessage("Database created OK");
}

void
CShimmerMainWindow::database_test(void)
{
    QFileDialog fileDialog;
    fileDialog.setFileMode(QFileDialog::AnyFile);
    fileDialog.setNameFilter(tr("SQLite3 DB (*.db)"));
    fileDialog.setDefaultSuffix(QString("db"));

    if (not fileDialog.exec())
        return;

    QStringList filenames = fileDialog.selectedFiles();
    if (filenames.size() == 0)
        return;

    std::string filename =
        filenames[ filenames.size()-1 ].toStdString();
    
    shimmer::infrastructure infra;

    int err = shimmer::load(filename, infra);
    if (err != SHIMMER_SUCCESS) {
        statusBar()->showMessage("Problem detected while loading DB");
        return;
    }

    statusBar()->showMessage("Database check OK");
}

void
CShimmerMainWindow::simulation_run(void)
{
    QFileDialog fileDialog;
    fileDialog.setFileMode(QFileDialog::AnyFile);
    fileDialog.setNameFilter(tr("SQLite3 DB (*.db)"));
    fileDialog.setDefaultSuffix(QString("db"));

    if (not fileDialog.exec())
        return;

    QStringList filenames = fileDialog.selectedFiles();
    if (filenames.size() == 0)
        return;

    std::string filename =
        filenames[ filenames.size()-1 ].toStdString();
    
    shimmer::config cfg;
    cfg.database = filename;
    launch_solver(cfg);

    statusBar()->showMessage("Simulation run OK");
}