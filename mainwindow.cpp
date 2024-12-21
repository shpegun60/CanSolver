#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <iostream>
#include "CanSolver.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    CanSolverTest();
}

MainWindow::~MainWindow()
{
    delete ui;
}
