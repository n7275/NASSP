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
#pragma once

#include <vector>
#include <string>
#include <fstream>
#include "Windows.h"

#define RTCC_TLM_START_STRING	"TLM_RTCC_BEGIN"
#define RTCC_TLM_END_STRING	    "TLM_RTCC_END"

class RTCC;
class MCC;
struct GroundStation;

namespace RTCC_Telemetry
{
	struct RTCC_TCP_TLM_Config
	{
		int CSM_PORT;
		int LEM_PORT; 
		int SIVbIU_PORT; 

		char* CSMDescriptorTable;
		char* CMCDescriptorTable;
		char* AMDDescriptorTable;

		char* LMDescriptorTable;
		char* LGCDescriptorTable;
		char* AGSDescriptorTable;

		char* SIVbIUDescriptorTable;
		char* SIIDescriptorTable;
		char* SICDescriptorTable;
	};

	class NASCOM
	{
	public:
		void refresh(double dt);
		inline void setActStation(GroundStation* station) { ActiveStation = station; };
	private:
		GroundStation* ActiveStation;
		double SignalStrength;
	};

	enum status
	{
		Missing,
		OffScaleHigh,
		OffscaleLow,
		Dead,
		Static,
		Nochange,
		Limit
	};

	enum TelemetryMeasurementTypes
	{
		RTCC_TLM_NULL,
		RTCC_TLM_A,
		RTCC_TLM_DP,
		RTCC_TLM_DS,
		RTCC_TLM_E,
		RTCC_TLM_SRC
	};

	enum TelemetryParameterUnits
	{
		UnitsNULL,
		Percentage,
		Sci,			// Revisit the validity of this one
		PSIA,
		PSIG,
		TempF,
		Volts,
		Amperes,
		PPH,
		DBM,
		degrees,
		G
	};

	struct ParameterAndStatus
	{
		char name[64];
		double parameter;
		bool parameterBit;
		status ParameterStatus;
	};

	struct DescriptorTableParameter
	{
		char name[64];
		unsigned int offset;
		TelemetryMeasurementTypes Type;
		unsigned int channel;
		unsigned int ccode;
		TelemetryParameterUnits Unit;
		double low;
		double high;
	};

	class DescriptorTableFormat
	{
	public:
		DescriptorTableFormat();
		void InitParameter(char* name, unsigned int offset, TelemetryMeasurementTypes Type, unsigned int channel, unsigned int ccode, TelemetryParameterUnits Unit, double low, double high); //might make this a private member later...
	private:
		std::vector<DescriptorTableParameter> Parameters;
	};

	class IntermediateDataArray
	{
	public:
		IntermediateDataArray(const unsigned int VEHCode, DescriptorTableFormat* Format);
		double GetParameter(char* name);
		double GetStatus(char* name);
	private:
		unsigned int vehicleIdentCode;
		GroundStation* TelemetrySite;
		bool live;
		bool highBit;
		double GMTA;
		double GMTR;

		std::vector<ParameterAndStatus> Data;
	};

	class TelemetryWorker
	{
	public:
		TelemetryWorker(const unsigned int WorkerSocket, const unsigned int VEHCode);
		~TelemetryWorker();

		int lock_type;	
		int frame_addr;
		int framect;
		int agc_lock_type;
		int agc_frame_addr;
		int agc_framect;

		IntermediateDataArray * intermediatedataarrays[3];

		// Winsock
		WSADATA wsaData;
		SOCKET m_socket;
		sockaddr_in clientService;
		int conn_status;

		void WinsockInit();
		void ConnectToHost();
		void CommThread();
		void parse_lbr(uint8_t recvdWord, int offset);
		void parse_hbr(uint8_t recvdWord, int offset);
	private:
		int SYNCWORDS[3];
		int LBRSYNC;
		int HBRSYNC;
		int LBRWORDCOUNT;
		int LBRFRAMECOUNT;
		int HBRWORDCOUNT;
		int HBRFRAMECOUNT;
	};

	class TelemetryProcessor
	{
	public:
		TelemetryProcessor();
		void Init(MCC *M, RTCC *R);
		void InitWorker(const unsigned int WorkerSocket, const char* DescriptorTableFile1, const char* DescriptorTableFile2, const char* DescriptorTableFile3);
		void ParseDescriptorTableFormatFile(char *file);
	private:
		MCC *mcc;
		RTCC *rtcc;
		std::vector<TelemetryWorker> Workers;
		std::vector<IntermediateDataArray> IntermediateDataArrays;
		std::vector<DescriptorTableFormat> DescriptorTableFormats;
		unsigned int sockets[64];
	};
}

bool papiReadTCPTelemetryConfigFile(char *line, char* item, RTCC_Telemetry::RTCC_TCP_TLM_Config &TCPtlmConfig);

/*
	Sources
	[1]	19680004357		PRELIMINARY PERFORMANCE ANALYSIS OF HIGH-SPEED DIGITAL DATA CIRCUITS IN THE NASCOM NETWORK
	[2] 19700024253		AS-508 MCC/MSFN MISSION CONFIGURATION SYSTEM DESCRIPTION
*/