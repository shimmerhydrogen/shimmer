#pragma once

#include <QtGui>
#include <QMainWindow>
#include <QMenu>
#include <iostream>

class CShimmerCentralWidget : public QWidget {
    Q_OBJECT

public:
    CShimmerCentralWidget(QWidget *parent);
    
};

class CShimmerMainWindow : public QMainWindow {
    Q_OBJECT

    QMenu *fileMenu_;
    QMenu *simMenu_;
    QMenu *dbMenu_;

private slots:
    void database_init(void);
    void database_test(void);
    void simulation_run(void);

public:
    CShimmerMainWindow(QWidget *parent = nullptr);

};