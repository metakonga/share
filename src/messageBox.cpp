#include "messageBox.h"

QMessageBox *messageBox::msg = NULL;

messageBox::messageBox()
{

}

messageBox::~messageBox()
{

}

int messageBox::run(QString text, QString info /* = "" */, QMessageBox::StandardButtons buttons /* = QMessageBox::Default */, QMessageBox::StandardButton button /* = QMessageBox::Default */)
{
	msg = new QMessageBox;
	msg->setText(text);
	if (info != "")
		msg->setInformativeText(info);
	if (buttons != QMessageBox::Default)
		msg->setStandardButtons(buttons);
	if (button != QMessageBox::Default)
		msg->setDefaultButton(button);
	int ret = msg->exec();
	delete msg;
	return ret;
}