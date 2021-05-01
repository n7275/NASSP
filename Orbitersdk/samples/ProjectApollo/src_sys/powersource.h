/***************************************************************************
  This file is part of Project Apollo - NASSP
  Copyright 2004-2005 Mark Grant

  ORBITER vessel module: Power distribution code

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

  See http://nassp.sourceforge.net/license/ for more details.

  **************************************************************************/

#if !defined(_PA_POWERSOURCE_H)
#define _PA_POWERSOURCE_H

#include "PanelSDK/PanelSDK.h"
#include "PanelSDK/Internals/ESystems.h"

namespace powerMergeCalc
{
	//This was generated using the following matlab script
	// clear;
	// clc;
	// pkg load symbolic;
	// sympref display unicode;

	// syms R0 R1 R2 R3 R4 R5 real;
	// syms V1 V2 V3 V4 V5 real;

	// ImpedanceMatrix5 = [R1+R0,  -R1,      0,      0,   0;...
	//                      -R1, R1+R2,    -R2,      0,   0;...
	//                        0,   -R2,  R2+R3,    -R3,   0;
	//                        0,     0,    -R3,  R3+R4,  -R4;
	//                        0,     0,      0,    -R4,  R4+R5];


	// Res5 = inv(ImpedanceMatrix5)*[V1,V2-V1,V3-V2,V4-V3,V5-V4]';

	// I1 = simplify(Res5(1));
	// I2 = simplify(Res5(2));
	// I3 = simplify(Res5(3));
	// I4 = simplify(Res5(4));
	// I5 = simplify(Res5(5));

	// P51 = ccode(simplify((I1-I2)*V1));
	// P52 = ccode(simplify((I2-I3)*V2));
	// P53 = ccode(simplify((I3-I4)*V3));
	// P54 = ccode(simplify((I4-I5)*V4));
	// P55 = ccode(simplify((I5-0)*V5));

	// clear I1 I2 I3 I4 I5;

	// ImpedanceMatrix4 = ImpedanceMatrix5(1:end-1,1:end-1);
	// Res4 = inv(ImpedanceMatrix4)*[V1,V2-V1,V3-V2,V4-V3]';
	// I1 = simplify(Res4(1));
	// I2 = simplify(Res4(2));
	// I3 = simplify(Res4(3));
	// I4 = simplify(Res4(4));

	// P41 = ccode(simplify((I1-I2)*V1));
	// P42 = ccode(simplify((I2-I3)*V2));
	// P43 = ccode(simplify((I3-I4)*V3));
	// P44 = ccode(simplify((I4-0)*V4));

	// clear I1 I2 I3 I4;

	// ImpedanceMatrix3 = ImpedanceMatrix4(1:end-1,1:end-1);
	// Res3 = inv(ImpedanceMatrix3)*[V1,V2-V1,V3-V2]';
	// I1 = simplify(Res3(1));
	// I2 = simplify(Res3(2));
	// I3 = simplify(Res3(3));

	// P31 = ccode(simplify((I1-I2)*V1));
	// P32 = ccode(simplify((I2-I3)*V2));
	// P33 = ccode(simplify((I3-0)*V3));

	// clear I1 I2 I3;

	// ImpedanceMatrix2 = ImpedanceMatrix3(1:end-1,1:end-1);
	// Res2 = inv(ImpedanceMatrix2)*[V1,V2-V1]';
	// I1 = simplify(Res2(1));
	// I2 = simplify(Res2(2));

	// P21 = ccode(simplify((I1-I2)*V1));
	// P22 = ccode(simplify(I2)*V2));

	inline void twoWay(double V1, double V2, double R0, double R1, double R2, double &P1, double &P2)
	{
		P1 = V1 * (R0*V1 - R0 * V2 + R2 * V1) / (R0*R1 + R0 * R2 + R1 * R2);
		P2 = V2 * (-R0 * V1 + R0 * V2 + R1 * V2) / (R0*R1 + R0 * R2 + R1 * R2);
	}

	inline void threeWay(double V1, double V2, double V3, double R0, double R1, double R2, double R3, double &P1, double &P2, double &P3)
	{
		double denominator = (R0*R1*R2 + R0 * R1*R3 + R0 * R2*R3 + R1 * R2*R3);
		P1 = V1 * (R0*R2*V1 - R0 * R2*V3 + R0 * R3*V1 - R0 * R3*V2 + R2 * R3*V1) / denominator;
		P2 = V2 * (R0*R1*V2 - R0 * R1*V3 - R0 * R3*V1 + R0 * R3*V2 + R1 * R3*V2) / denominator;
		P3 = V3 * (-R0 * R1*V2 + R0 * R1*V3 - R0 * R2*V1 + R0 * R2*V3 + R1 * R2*V3) / denominator;
	}

	inline void fourWay(double V1, double V2, double V3, double V4, double R0, double R1, double R2, double R3, double R4, double &P1, double &P2, double &P3, double &P4)
	{
		double denominator = (R0*R1*R2*R3 + R0 * R1*R2*R4 + R0 * R1*R3*R4 + R0 * R2*R3*R4 + R1 * R2*R3*R4);
		P1 = V1 * (R0*R2*R3*V1 - R0 * R2*R3*V4 + R0 * R2*R4*V1 - R0 * R2*R4*V3 + R0 * R3*R4*V1 - R0 * R3*R4*V2 + R2 * R3*R4*V1) / denominator;
		P2 = V2 * (R0*R1*R3*V2 - R0 * R1*R3*V4 + R0 * R1*R4*V2 - R0 * R1*R4*V3 - R0 * R3*R4*V1 + R0 * R3*R4*V2 + R1 * R3*R4*V2) / denominator;
		P3 = V3 * (R0*R1*R2*V3 - R0 * R1*R2*V4 - R0 * R1*R4*V2 + R0 * R1*R4*V3 - R0 * R2*R4*V1 + R0 * R2*R4*V3 + R1 * R2*R4*V3) / denominator;
		P4 = V4 * (-R0 * R1*R2*V3 + R0 * R1*R2*V4 - R0 * R1*R3*V2 + R0 * R1*R3*V4 - R0 * R2*R3*V1 + R0 * R2*R3*V4 + R1 * R2*R3*V4) / denominator;
	}

	inline void fiveWay(double V1, double V2, double V3, double V4, double V5, double R0, double R1, double R2, double R3, double R4, double R5, double &P1, double &P2, double &P3, double &P4, double &P5)
	{
		double denominator = (R0*R1*R2*R3*R4 + R0 * R1*R2*R3*R5 + R0 * R1*R2*R4*R5 + R0 * R1*R3*R4*R5 + R0 * R2*R3*R4*R5 + R1 * R2*R3*R4*R5);
		P1 = V1 * (R0*R2*R3*R4*V1 - R0 * R2*R3*R4*V5 + R0 * R2*R3*R5*V1 - R0 * R2*R3*R5*V4 + R0 * R2*R4*R5*V1 - R0 * R2*R4*R5*V3 + R0 * R3*R4*R5*V1 - R0 * R3*R4*R5*V2 + R2 * R3*R4*R5*V1) / denominator;
		P2 = V2 * (R0*R1*R3*R4*V2 - R0 * R1*R3*R4*V5 + R0 * R1*R3*R5*V2 - R0 * R1*R3*R5*V4 + R0 * R1*R4*R5*V2 - R0 * R1*R4*R5*V3 - R0 * R3*R4*R5*V1 + R0 * R3*R4*R5*V2 + R1 * R3*R4*R5*V2) / denominator;
		P3 = V3 * (R0*R1*R2*R4*V3 - R0 * R1*R2*R4*V5 + R0 * R1*R2*R5*V3 - R0 * R1*R2*R5*V4 - R0 * R1*R4*R5*V2 + R0 * R1*R4*R5*V3 - R0 * R2*R4*R5*V1 + R0 * R2*R4*R5*V3 + R1 * R2*R4*R5*V3) / denominator;
		P4 = V4 * (R0*R1*R2*R3*V4 - R0 * R1*R2*R3*V5 - R0 * R1*R2*R5*V3 + R0 * R1*R2*R5*V4 - R0 * R1*R3*R5*V2 + R0 * R1*R3*R5*V4 - R0 * R2*R3*R5*V1 + R0 * R2*R3*R5*V4 + R1 * R2*R3*R5*V4) / denominator;
		P5 = V5 * (-R0 * R1*R2*R3*V4 + R0 * R1*R2*R3*V5 - R0 * R1*R2*R4*V3 + R0 * R1*R2*R4*V5 - R0 * R1*R3*R4*V2 + R0 * R1*R3*R4*V5 - R0 * R2*R3*R4*V1 + R0 * R2*R3*R4*V5 + R1 * R2*R3*R4*V5) / denominator;
	}
}

class PowerSource : public e_object {

public:
	PowerSource();
	~PowerSource();

	double Voltage();
	double Current();

protected:
};

///
/// \ingroup PowerSource
/// \brief Allows drawing power from two sources simultaneously. Optionally simulates diodes on none, either, or both of the inputs.
/// 

class PowerMerge : public PowerSource{
public:
	PowerMerge(char *i_name, PanelSDK &p, bool HasDiodeA, bool HasDiodeB, double ISatA, double ISatB, double InResistenceA, double InResistenceB);
	PowerMerge(char *i_name, PanelSDK &p);
	
	//Methods derived From e_object
	void DrawPower(double watts);
	double Voltage();
	double Current();
	void UpdateFlow(double dt);

	//PowerMerge methods
	void WireToBuses(e_object *a, e_object *b);
	

protected:
	PanelSDK &sdk;

	e_object *BusA;
	e_object *BusB;
	bool IsDiodeA, IsDiodeB;
	double R1, R2;
	double VoltsA, VoltsB;
	double Power1, Power2;
	Diode DiodeA, DiodeB;
	e_object *DrawA;
	e_object *DrawB;
};

///
/// \ingroup PowerSource
/// \brief Allows drawing power from three sources simultaneously. Optionally simulates diodes on any combination of inputs inputs.
///	\see PowerMerge
/// 

class ThreeWayPowerMerge : public PowerSource {
public:
	ThreeWayPowerMerge(char *i_name, PanelSDK &p);
	double Voltage();
	void DrawPower(double watts);
	void WireToBuses(e_object *a, e_object *b, e_object *c) { Phase1 = a; Phase2 = b; Phase3 = c; };
	void WireToBus(int bus, e_object* e);
	bool IsBusConnected(int bus);
	double Current();

protected:
	PanelSDK &sdk;

	e_object *Phase1;
	e_object *Phase2;
	e_object *Phase3;
};

class NWayPowerMerge : public PowerSource {
public:
	NWayPowerMerge(char *i_name, PanelSDK &p, int nSources);
	~NWayPowerMerge();

	double Voltage();
	void DrawPower(double watts);
	void WireToBus(int bus, e_object* e);
	bool IsBusConnected(int bus);
	double Current();

protected:
	PanelSDK &sdk;

	int nSources;
	e_object **sources;
};

class PowerBreaker : public PowerSource {

public:
	PowerBreaker() { breaker_open = false; };
	double Voltage();
	bool IsOpen() { return breaker_open; };
	virtual void SetOpen(bool state) { breaker_open = state; };

protected:
	bool breaker_open;
};

class PowerSDKObject : public PowerSource {

public:
	PowerSDKObject() { SDKObj = 0; };

	double Voltage();
	double Current();
	double PowerLoad();

	void DrawPower(double watts);
	void WireToSDK(e_object *s) { SDKObj = s; };

protected:
	e_object *SDKObj;
};


class DCBusController : public e_object {

public:
	DCBusController(char *i_name, PanelSDK &p);
	void Init(e_object *fc1, e_object *fc2, e_object *fc3, e_object *bat1, e_object *bat2, e_object *gse, e_object *bc1, e_object *bc2, e_object *bc3);
	void refresh(double dt);
	void ConnectFuelCell(int fc, bool connect);
	bool IsFuelCellConnected(int fc);
	bool IsBusContPowered(int fc);
	bool IsSMBusPowered();
	bool IsFuelCellDisconnectAlarm();
	e_object *GetBusSource() { return &busPower; };
	void SetTieState(int s) { tieState = s; };
	void SetTieAuto(bool s) { tieAuto = s; };
	void SetGSEState(int s);
	void Load(char *line);
	void Save(FILEHANDLE scn);

protected:
	PanelSDK &sdk;

	e_object *fuelcell1, *fuelcell2, *fuelcell3;
	e_object *battery1, *battery2;
	e_object *gseBattery;
	e_object *busCont1, *busCont2, *busCont3;

	NWayPowerMerge busPower;
	ThreeWayPowerMerge fcPower;
	PowerMerge batPower;
	bool fcDisconnectAlarm[4];
	int tieState;
	bool tieAuto;
	int gseState;
};


class BatteryCharger : public e_object {

public:
	BatteryCharger(char *i_name, PanelSDK &p);
	void Init(e_object *bat1, e_object *bat2, e_object *bat3, 
		      e_object *batSup1, e_object *batSup2, e_object *batSup3,
			  e_object *dc1, e_object *dc2, e_object* ac,
			  e_object *bat3pwr);
	void UpdateFlow(double dt);
	void DrawPower(double watts);
	void Charge(int bat);
	void Load(char *line);
	void Save(FILEHANDLE scn);

protected:
	PanelSDK &sdk;

	e_object *battery1, *battery2, *battery3;
	e_object *batSupply1, *batSupply2, *batSupply3;

	PowerMerge dcPower;
	e_object *acPower;

	e_object *currentBattery;

	e_object *bat3Power;
};

class Connector;
///
/// \ingroup PowerSource
/// \brief Connector object: sends messages to a connector, rather than to the next e_object
/// in the chain. Mostly used for connecting to another vessel.
///
class PowerSourceConnectorObject : public e_object
{
public:
	PowerSourceConnectorObject(char *i_name, PanelSDK &p);
	void SetConnector(Connector *c) { connect = c; };

	double Voltage();
	double Current();
	void UpdateFlow(double dt);
	void DrawPower(double watts);

protected:
	Connector *connect;
};

///
/// \brief Connector object: receives messages from a connector, rather than the next e_object
/// in the chain. Mostly used for connecting to another vessel.
///
class PowerDrainConnectorObject : public e_object
{
public:
	PowerDrainConnectorObject(char *i_name, PanelSDK &sdk);

	void SetConnector(Connector *c) { connect = c; };

	void ProcessUpdateFlow(double dt);
	void ProcessDrawPower(double watts);
	void refresh(double dt);
	void Disconnected();

protected:
	Connector *connect;

	double PowerDraw;
	double PowerDrawn;
};

///
/// \ingroup Connectors
/// \brief Message type to process electric power through a connector.
///
enum PowerSourceMessageType
{
	POWERCON_GET_VOLTAGE,					///< Get voltage.
	POWERCON_GET_CURRENT,					///< Get current.
	POWERCON_DRAW_POWER,					///< Draw power from connector.
	POWERCON_UPDATE_FLOW,					///< Update power flow.
};

#endif
