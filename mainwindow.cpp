#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <iostream>
#include "CanSolver.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    CanSolver nnn;
    CanSolver::Input nominal;
    CanSolver::calculate4_can(nominal);
    CanSolver::test();
}

MainWindow::~MainWindow()
{
    delete ui;
}
