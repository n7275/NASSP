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

   void RTCC_Telemetry::DownlistFormat::InitParameter(char* name, unsigned int offset, TelemetryMeasurementTypes Type, unsigned int channel, unsigned int ccode, TelemetryParameterUnits Unit, double low, double high)
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
	int bytesRecv = SOCKET_ERROR;
	uint8_t recvbuf[1024];
	int die = 0;

	int frame_count = 0;
	int word_addr = 0;
	int byte_offset = 0;
	int bytect = 0;

	lock_type = 0;   agc_lock_type = 0;
	frame_addr = 0;	 agc_frame_addr = 0;
	framect = 0;     agc_framect = 0;
	
	while (!die)
	{
		bytesRecv = recv(m_socket, (char *)recvbuf, 1024, 0);
		if (bytesRecv == SOCKET_ERROR)
		{
			closesocket(m_socket);
			conn_status = 0;
			return;
		}
		else
		{
			byte_offset = 0;
			while (byte_offset < bytesRecv)
			{
				switch (lock_type)
				{
					case 0: // OUT SYNC 0
						bytect = 0;
						if (recvbuf[byte_offset] == SYNCWORDS[0])  // Sync char 1 recieved
						{
							lock_type = 1;
						}
						break;

					case 1: // OUT SYNC 1
						if (recvbuf[byte_offset] == SYNCWORDS[1]) // Sync char 2 recieved
						{
							lock_type = 2;
						}
						else
						{
							lock_type = 0;
						}
						break;

					case 2: // OUT SYNC 2
						if (recvbuf[byte_offset] == SYNCWORDS[2]) // Sync char 3 recieved
						{
							lock_type = 3;
						}
						else
						{
							lock_type = 0;
						}
						break;
					case 3: // OUT SYNC 3
						framect = recvbuf[byte_offset] & 077;
						lock_type++;
						break;
					case 4: // OUT SYNC 4 & CHECK FOR HBR OR LBR
						if (recvbuf[byte_offset] == LBRSYNC) // LBR Sync char 4 recieved
						{
							lock_type = 10;
							bytect = 5;
						}
						else if (recvbuf[byte_offset] == HBRSYNC)
						{
							lock_type = 20;
							bytect = 5;
						}
						else
						{
							lock_type = 0;
						}
						break;
					case 10: //LBR SYNC
						parse_lbr(recvbuf[byte_offset], bytect);
						bytect++;
						if (bytect > LBRWORDCOUNT-1)
						{
							bytect = 0;
						}
						break;
					case 20: // HBR
						parse_hbr(recvbuf[byte_offset], bytect);
						bytect++;
						if (bytect > HBRWORDCOUNT-1)
						{
							bytect = 0;
							frame_addr++;
							if (frame_addr > HBRFRAMECOUNT-1)
							{
								frame_addr = 0;
							}
						}
						break;
				}
				byte_offset++;
			}
		}
	}
}

void RTCC_Telemetry::TelemetryWorker::parse_lbr(uint8_t recvdWord, int offset)
{
}

void RTCC_Telemetry::TelemetryWorker::parse_hbr(uint8_t recvdWord, int offset)
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
