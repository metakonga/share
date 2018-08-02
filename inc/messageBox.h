#ifndef MESSAGEBOX_H
#define MESSAGEBOX_H

#include <QMessageBox>
#include <QString>

class messageBox
{
public:
	messageBox();
	~messageBox();

	static int run(QString text, QString info = "", QMessageBox::StandardButtons buttons = QMessageBox::Default, QMessageBox::StandardButton button = QMessageBox::Default);
private:
	static QMessageBox *msg;
};

#endif