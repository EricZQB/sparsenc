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
extern "C" {
    #include "sparsenc.h"   // libsparsenc is written in C
}

//
//
//
// Build network topology
// =====================================================================
//
//                    _ _ _ _ _ _ _ _ _ _ _ _ _ _
//         p-2-p     |             |             |      p-2-p
//   n0 ------------ | NetDevice1  |  NetDevice2 | ------------- n2
//                   |_ _ _ _ _ _ _|_ _ _ _ _ _ _|               
//                                 n1
//
// (Encoder)                    (Recoder)                     (Decoder)
//
// ======================================================================

using namespace ns3;

typedef struct snc_parameters   SncParameter;
typedef struct snc_packet       SncPacket;
typedef struct snc_context      SncEncoder;
typedef struct snc_buffer_bats  SncRecoder;
typedef struct snc_decoder      SncDecoder;

NS_LOG_COMPONENT_DEFINE ("SparseNcApplication");

class SparseNcApplication
{
public:
  SparseNcApplication(SncParameter sncParameter);
  void EncoderSend(Ptr<Socket> socket);
  void DecoderRecv(Ptr<Socket> socket);
  void RecoderRecv(Ptr<Socket> socket);
  void RecoderSend(Ptr<Socket> socket);
  void SetEncoderRate(DataRate dataRate);
  void SetRecoderRate(DataRate dataRate);

private:
  SncParameter m_sp;
  SncEncoder  *m_encoder;
  SncDecoder  *m_decoder;
  SncRecoder  *m_recoder;
  DataRate     m_encoderRate;
  DataRate     m_recoderRate;
};

SparseNcApplication::SparseNcApplication(SncParameter sncParameter)
{
  srand(static_cast<uint32_t>(time(0)));
  m_sp = sncParameter;
  // A buffer containing random bytes
  unsigned char *buf = (unsigned char *) malloc(m_sp.datasize);
  for (int i=0; i<m_sp.datasize; i++) {
    buf[i] = rand() % (1<<8);
  }
  m_encoder = snc_create_enc_context(buf, &m_sp);
  m_sp.seed = (snc_get_parameters(m_encoder))->seed;
  m_decoder = snc_create_decoder(&m_sp, CBD_DECODER);
  uint32_t bufferSize = 10;
  m_recoder = snc_create_buffer_bats(&m_sp, bufferSize);
}

void 
SparseNcApplication::EncoderSend(Ptr<Socket> socket)
{
  if (!snc_decoder_finished(m_decoder))
    {
      SncPacket *ncpkt = snc_generate_packet(m_encoder);
      unsigned char *pktstr = snc_serialize_packet(ncpkt, &m_sp);
      snc_free_packet(ncpkt);
      int nbytes = snc_packet_length(&m_sp);
      // std::cout << "SNC packet size: " << nbytes << std::endl;
      Ptr<Packet> packet = Create<Packet>(pktstr, nbytes);
      socket->Send(packet);
      NS_LOG_UNCOND (
        "[Encoder] At time " 
        << Simulator::Now ().GetSeconds () 
        << " (s) server sent a packet of size "
        << packet->GetSize()
        << " bytes"
      );
      // schedule the next transmission
      Time tNext (Seconds (snc_packet_length(&m_sp) * 8 / static_cast<double> (m_encoderRate.GetBitRate ())));
      Simulator::Schedule(tNext, &SparseNcApplication::EncoderSend, this, socket);
    } 
  else 
    {
      NS_LOG_UNCOND (
        "At time "
        << Simulator::Now ().GetSeconds () 
        << " (s) client completes decoding! "
      );
      print_code_summary(m_encoder, snc_decode_overhead(m_decoder), snc_decode_cost(m_decoder));
      // verify decoded data
      unsigned char *data_o = snc_recover_data(m_encoder);
      unsigned char *data_d = snc_recover_data(snc_get_enc_context(m_decoder));
      if (memcmp(data_o, data_d, m_sp.datasize) != 0) {
        std::cout << "ERROR: recovered is NOT identical to original" << std::endl;
      } else {
        std::cout << "All source packets are recovered correctly" << std::endl;
      }
      socket->Close ();
    }
}

void 
SparseNcApplication::DecoderRecv(Ptr<Socket> socket)
{
    Ptr<Packet> pkt = socket->Recv();   // receive serialized snc_packet from socket
    unsigned char *pktstr = (unsigned char *) calloc(snc_packet_length(&m_sp), sizeof(unsigned char));
    pkt->CopyData(pktstr, snc_packet_length(&m_sp));
    SncPacket *ncpkt = snc_deserialize_packet(pktstr, &m_sp);
    snc_process_packet(m_decoder, ncpkt);
    NS_LOG_UNCOND (
      "[Decoder] At time "
      << Simulator::Now ().GetSeconds () 
      << " (s) client received a packet "
    );
}


void 
SparseNcApplication::RecoderRecv(Ptr<Socket> socket)
{
    Ptr<Packet> packet = socket->Recv();   // receive serialized snc_packet from socket
    unsigned char *pktstr = (unsigned char *) calloc(snc_packet_length(&m_sp), sizeof(unsigned char));
    packet->CopyData(pktstr, snc_packet_length(&m_sp));
    SncPacket *sncPacket = snc_deserialize_packet(pktstr, &m_sp);
    snc_buffer_packet_bats(m_recoder, sncPacket);   // intermediate ndoes buffer packets
    NS_LOG_UNCOND (
      "[Recoder] At time "
      << Simulator::Now ().GetSeconds () 
      << " (s) recoder buffered a packet "
    );
}


void 
SparseNcApplication::RecoderSend(Ptr<Socket> socket)
{
  if (!snc_decoder_finished(m_decoder)) 
    {
      SncPacket *sncPacket = snc_recode_packet_bats(m_recoder);
      if (sncPacket != NULL) {
        unsigned char *pktstr = snc_serialize_packet(sncPacket, &m_sp);
        snc_free_packet(sncPacket);
        Ptr<Packet> packet = Create<Packet>(pktstr, snc_packet_length(&m_sp));
        socket->Send(packet);
        NS_LOG_UNCOND (
          "[Recoder] At time "
          << Simulator::Now ().GetSeconds () 
          << " (s) recoder sent a recoded packet "
        );
      } else {
        NS_LOG_UNCOND (
          "[Recoder] At time "
          << Simulator::Now ().GetSeconds () 
          << " (s) recoder forgone a transmission opportunity "
        );
      }
      // schedule the next transmission
      Time tNext (Seconds (snc_packet_length(&m_sp) * 8 / static_cast<double> (m_recoderRate.GetBitRate ())));
      Simulator::Schedule(tNext, &SparseNcApplication::RecoderSend, this, socket);
    }
}

void 
SparseNcApplication::SetEncoderRate(DataRate dataRate)
{
  m_encoderRate = dataRate;
}

void 
SparseNcApplication::SetRecoderRate(DataRate dataRate)
{
  m_recoderRate = dataRate;
}

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
MacTxDrop (std::string context, Ptr<const Packet> p)
{
  NS_LOG_UNCOND ("MacTxDrop " << p->GetSize() << " bytes at " << Simulator::Now ().GetSeconds () << " << "<< context);
}

static void
UdpDrop (std::string context, Ptr<const Packet> p)
{
  NS_LOG_UNCOND ("UdpDrop " << p->GetSize() << " bytes at " << Simulator::Now ().GetSeconds () << " << "<< context);
}

int 
main (int argc, char *argv[])
{

  SncParameter  sp;
  sp.type     = BATS_SNC;
  sp.datasize = 168000;
  sp.size_p   = 1050;
  sp.bpc      = 1;
  sp.sys      = 0;
  sp.seed     = -1;

  uint32_t gfPower   = 8;
  uint32_t batchSize = 160;
  uint32_t bts       = 1000;
  double   lossRate  = 0.01;
  CommandLine cmd;
  cmd.AddValue ("gfPower",
                "Galois field size as a power of 2 (default gfPower=8, i.e., q=2^8)", gfPower);
  cmd.AddValue ("batchSize",
                "Batch size of FD-BATS code", batchSize);
  cmd.AddValue ("bts",
                "Batch transmission size of FD-BATS code", bts);
  cmd.AddValue ("lossRate",
                "Packet loss rate over the link (loss rate identical on two hops) ", lossRate);
  cmd.Parse (argc, argv);

  Time::SetResolution (Time::NS);
  
  sp.size_c   = batchSize == 160 ? 0 : 23;  // number of precoding packets
  sp.size_g   = batchSize;
  sp.size_b   = bts;
  sp.gfpower  = gfPower;

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

  SparseNcApplication app = SparseNcApplication (sp);      // A network coding application instance
  app.SetEncoderRate (DataRate ("10Mbps"));
  app.SetRecoderRate (DataRate ("10Mbps"));
  
  uint32_t port = 6666;
  TypeId tid = TypeId::LookupByName ("ns3::UdpSocketFactory");
  // Encoder: as a UDP server on node 0
  Ptr<Socket> encoderSocket = Socket::CreateSocket(firstHopNodes.Get(0), tid);
  InetSocketAddress recoderAddr = InetSocketAddress (firstHopInterfaces.GetAddress (1), port);
  encoderSocket->Connect (recoderAddr);

  // Recoder part 1: as a UDP client for receiving from encoder on node 1
  Ptr<Socket> recoderSocketIn = Socket::CreateSocket(firstHopNodes.Get(1), tid);
  InetSocketAddress recoderLocal = InetSocketAddress (Ipv4Address::GetAny (), port);
  recoderSocketIn->Bind (recoderLocal);
  recoderSocketIn->SetRecvCallback (MakeCallback (&SparseNcApplication::RecoderRecv, &app));
  
  // Recoder part 2: as a UDP server for sending recoded packet to decoder on node 1
  Ptr<Socket> recoderSocketOut = Socket::CreateSocket(secondHopNodes.Get(0), tid);
  InetSocketAddress decoderAddr = InetSocketAddress (secondHopInterfaces.GetAddress (1), port);
  recoderSocketOut->Connect (decoderAddr);

  // Decoder: snc decoder as a UDP client (i.e., receiving packets) on node 2
  Ptr<Socket> decoderSocket = Socket::CreateSocket (secondHopNodes.Get (1), tid);
  InetSocketAddress decoderLocal = InetSocketAddress (Ipv4Address::GetAny (), port);
  decoderSocket->Bind (decoderLocal);
  decoderSocket->SetRecvCallback (MakeCallback (&SparseNcApplication::DecoderRecv, &app));

  Config::Connect ("/NodeList/*/DeviceList/*/$ns3::PointToPointNetDevice/MacRx", MakeCallback (&MacRx));
  Config::Connect ("/NodeList/*/DeviceList/*/$ns3::PointToPointNetDevice/MacTx", MakeCallback (&MacTx));
  Config::Connect ("/NodeList/*/DeviceList/*/$ns3::PointToPointNetDevice/MacTxDrop", MakeCallback (&MacTxDrop));
  Config::Connect ("/NodeList/*/DeviceList/*/$ns3::UdpSocketImpl/Drop", MakeCallback (&UdpDrop));

  // Schedule the first encoded packet transmission
  Simulator::ScheduleWithContext (encoderSocket->GetNode ()->GetId (), 
                                  Seconds (0.0), 
                                  &SparseNcApplication::EncoderSend, 
                                  &app,
                                  encoderSocket);

  // Schedule the first recoded packet transmission
  Simulator::ScheduleWithContext (recoderSocketOut->GetNode ()->GetId (), 
                                  Seconds (0.0), 
                                  &SparseNcApplication::RecoderSend, 
                                  &app,
                                  recoderSocketOut);

  Simulator::Stop (Seconds (100.0));
  Simulator::Run();
  Simulator::Destroy();

  return 0;
}