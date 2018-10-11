#include "errors.h"
#include "messageBox.h"

errors::error_type errors::e = errors::NO_ERROR;

errors::errors()
{

}

errors::~errors()
{

}

void errors::Error(QString &target)
{
	QString msg;
	QTextStream qts(&msg);
	switch (e)
	{
	case MBD_EXCEED_NR_ITERATION:
		qts << "Fail model : " + target << endl
			<< "Problem : Exceed the N-R iteration number." << endl
			<< "Tip : Check your multi-body model" << endl;
		break;
	}
	if (msg.isEmpty())
		return;
	messageBox::run(msg);
}