/***************************************************************************
  This file is part of Project Apollo - NASSP
  Copyright 2021

  MCC/RTCC Telemetry Classes

  Project Apollo is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  Project Apollo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with Project Apollo; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

  See https://github.com/orbiternassp/NASSP/blob/Orbiter2016/COPYING.txt
  for more details.

  **************************************************************************/
#include "RTCC_Telemetry.h"

void RTCC_Telemetry::DownlistFormat::InitParameter(char * name, TelemetryMeasurementTypes Type, unsigned int channel, unsigned int ccode, TelemetryParameterUnits Unit, double low, double high)
{
	RTCC_Telemetry::DownlistParameter TempParameter;

	strcpy(TempParameter.name, name);
	TempParameter.Type = Type;
	TempParameter.channel = channel;
	TempParameter.ccode = ccode;
	TempParameter.Unit = Unit;
	TempParameter.low = low;
	TempParameter.high = high;

	this->Parameters.push_back(TempParameter);
}

void RTCC_Telemetry::TelemetryProcessor::WinsockInit()
{
}

void RTCC_Telemetry::TelemetryProcessor::ConnectToHost()
{
}


void RTCC_Telemetry::TelemetryWorker::InitWorker()
{
}

void RTCC_Telemetry::TelemetryWorker::CommThread()
{
}

RTCC_Telemetry::IntermediateDataArray::IntermediateDataArray(unsigned int VEHCode, DownlistFormat * Format)
{
}

double RTCC_Telemetry::IntermediateDataArray::GetParameter(char * name)
{
	return 0.0;
}

double RTCC_Telemetry::IntermediateDataArray::GetStatus(char * name)
{
	return 0.0;
}
