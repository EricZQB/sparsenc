/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as
 * published by the Free Software Foundation;
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/point-to-point-module.h"
#include "ns3/applications-module.h"
#include "ns3/ipv4-global-routing-helper.h"

// Build network topology
// =====================================================================
//
//                    _ _ _ _ _ _ _ _ _ _ _ _ _ _
//         p-2-p     |             |             |      p-2-p
//   n0 ------------ | NetDevice1  |  NetDevice2 | ------------- n2
//                   |_ _ _ _ _ _ _|_ _ _ _ _ _ _|               
//                                 n1
//
// ======================================================================

using namespace ns3;

NS_LOG_COMPONENT_DEFINE ("LossyTwoHopTcp");

static void
RxDrop (Ptr<const Packet> p)
{
  NS_LOG_UNCOND ("RxDrop at " << Simulator::Now ().GetSeconds ());
}

static void
MacRx (std::string context, Ptr<const Packet> p)
{
  NS_LOG_UNCOND ("MacRx " << p->GetSize() << " bytes at " << Simulator::Now ().GetSeconds () << " << " << context);
}

static void
MacTx (std::string context, Ptr<const Packet> p)
{
  NS_LOG_UNCOND ("MacTx " << p->GetSize() << " bytes at " << Simulator::Now ().GetSeconds () << " << "<< context);
}


static void
BulkTx (std::string context, Ptr<const Packet> p)
{
  NS_LOG_UNCOND ("BulkTx " << p->GetSize() << " bytes at " << Simulator::Now ().GetSeconds ());
}

static void
SinkRx (std::string context, Ptr<const Packet> p, const Address &address)
{
  static int totalRxBytes = 0;
  totalRxBytes += p->GetSize();
  NS_LOG_UNCOND ("SinkRx a total of " << totalRxBytes << " bytes at " << Simulator::Now ().GetSeconds ());
}

/* main function */
int 
main (int argc, char *argv[])
{
  uint32_t maxBytes = 168000;
  double   lossRate = 0.01;
  CommandLine cmd;
  cmd.AddValue ("maxBytes",
                "Total number of bytes for application to send", maxBytes);
  cmd.AddValue ("lossRate",
                "Packet loss rate over the link (loss rate identical on two hops) ", lossRate);
  cmd.Parse (argc, argv);
  std::cout << "Loss rate: " << lossRate << std::endl;
  
  Time::SetResolution (Time::NS);

  // Create two nodes connected by the PPP link-layer protocol
  NodeContainer firstHopNodes;
  firstHopNodes.Create (2);

  PointToPointHelper pointToPoint;
  pointToPoint.SetDeviceAttribute ("DataRate", StringValue ("5Mbps"));
  pointToPoint.SetChannelAttribute ("Delay", StringValue ("50ms"));

  NetDeviceContainer firstHopDevices;
  firstHopDevices = pointToPoint.Install (firstHopNodes);

  // Add packet loss model and trace to the first hop
  Ptr<RateErrorModel> firstHopEm = CreateObject<RateErrorModel> ();
  firstHopEm->SetAttribute ("ErrorUnit", EnumValue (RateErrorModel::ERROR_UNIT_PACKET));
  firstHopEm->SetAttribute ("ErrorRate", DoubleValue (lossRate));
  firstHopDevices.Get (1)->SetAttribute ("ReceiveErrorModel", PointerValue (firstHopEm));
  firstHopDevices.Get (1)->TraceConnectWithoutContext ("PhyRxDrop", MakeCallback (&RxDrop));

  // Create two nodes for the second hop
  NodeContainer secondHopNodes;
  secondHopNodes.Add (firstHopNodes.Get (1));
  secondHopNodes.Create (1);

  // pointToPoint.SetDeviceAttribute ("DataRate", StringValue ("5Mbps"));
  // pointToPoint.SetChannelAttribute ("Delay", StringValue ("20ms"));

  NetDeviceContainer secondHopDevices;
  secondHopDevices = pointToPoint.Install (secondHopNodes);

  // Add packet loss model and trace to the second hop
  Ptr<RateErrorModel> secondHopEm = CreateObject<RateErrorModel> ();
  secondHopEm->SetAttribute ("ErrorUnit", EnumValue (RateErrorModel::ERROR_UNIT_PACKET));
  secondHopEm->SetAttribute ("ErrorRate", DoubleValue (lossRate));
  secondHopDevices.Get (1)->SetAttribute ("ReceiveErrorModel", PointerValue (secondHopEm));
  secondHopDevices.Get (1)->TraceConnectWithoutContext ("PhyRxDrop", MakeCallback (&RxDrop));

  InternetStackHelper stack;
  stack.Install (firstHopNodes.Get (0));
  stack.Install (secondHopNodes);

  // Assign IP addresses to the net devices
  Ipv4AddressHelper ipv4;
  // first hop is in 10.1.1.0
  ipv4.SetBase ("10.1.1.0", "255.255.255.0");
  Ipv4InterfaceContainer firstHopInterfaces = ipv4.Assign (firstHopDevices);
  // second hop is in 10.1.2.0
  ipv4.SetBase ("10.1.2.0", "255.255.255.0");
  Ipv4InterfaceContainer secondHopInterfaces = ipv4.Assign (secondHopDevices);

  uint32_t port = 6666;
  // Setup application (source) on node 0, i.e., firstHopNodes.Get (0)
  BulkSendHelper source ("ns3::TcpSocketFactory",
                         InetSocketAddress (secondHopInterfaces.GetAddress (1), port));
  // Set the amount of data to send in bytes.  Zero is unlimited.
  source.SetAttribute ("MaxBytes", UintegerValue (maxBytes));
  ApplicationContainer sourceApps = source.Install (firstHopNodes.Get (0));
  sourceApps.Start (Seconds (0.0));
  sourceApps.Stop (Seconds (100.0));
  //Config::Connect ("/NodeList/0/ApplicationList/0/$ns3::BulkSendApplication/Tx", MakeCallback (&BulkTx));

  // Setup application (sink) on node 2, i.e., secondHopNodes.Get (1)
  PacketSinkHelper sink ("ns3::TcpSocketFactory",
                         InetSocketAddress (Ipv4Address::GetAny (), port));
  ApplicationContainer sinkApps = sink.Install (secondHopNodes.Get (1));
  sinkApps.Start (Seconds (0.0));
  sinkApps.Stop (Seconds (100.0));
  Config::Connect ("/NodeList/2/ApplicationList/0/$ns3::PacketSink/Rx", MakeCallback (&SinkRx));

  Ipv4GlobalRoutingHelper::PopulateRoutingTables ();

  // Add PPP device trace
  Config::Connect ("/NodeList/*/DeviceList/*/$ns3::PointToPointNetDevice/MacRx", MakeCallback (&MacRx));
  Config::Connect ("/NodeList/*/DeviceList/*/$ns3::PointToPointNetDevice/MacTx", MakeCallback (&MacTx));

  NS_LOG_INFO ("Run Simulation.");
  Simulator::Stop (Seconds (100.0));
  Simulator::Run();
  Simulator::Destroy();

  Ptr<PacketSink> sink1 = DynamicCast<PacketSink> (sinkApps.Get (0));
  std::cout << "Total Bytes Received: " << sink1->GetTotalRx () << std::endl;

  return 0;
}