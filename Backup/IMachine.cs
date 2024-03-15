using System;
namespace Electroimpact
{
  namespace Machine
  {
    public interface IMachine
    {
      double A { get; }
      bool AddConfigFile(string FileName);
      double B { get; }
      double[] C { get; }
      bool CompTableOverride { get; set; }
      string[] CompTableVariableNames { get; }
      bool ContainsAxisAttribute(string AxName, string AttributeName);
      string[] DHStrings { get; }
      double[,] DotProduct();
      double[,] DotProduct(uint start, uint finish);
      string[] GetAttributeNames();
      string[] GetAxisNames();
      string[] GetAxisNames(bool ShowHidden);
      double GetAxisPostion(string axName);
      double GetAxisScaledPosition(string axName);
      System.Collections.Generic.List<Electroimpact.Machine.cCompTable.cCompStation> GetCompensationTable(string AxisName, double min, double max, out int lower, out int upper);
      System.Collections.Generic.List<Electroimpact.Machine.cCompTable.cCompStation> GetCompensationTable(string AxisName);
      string[] GetCompensationTableAxisNames { get; }
      int GetCompensationTableStationCount(string AxisName, double min, double max);
      int GetCompensationTableStationCount(string AxisName);
      bool IsHookedUp { get; }
      void NullCompTables();
      void OverrideCompTableValues(string CompTableVaiableName, double value);
      void PutCompensationTable(string AxisName, System.Collections.Generic.List<Electroimpact.Machine.cCompTable.cCompStation> Table);
      double ReadAttribute(string AttribName);
      double ReadAxisAttribute(string axName, string attributeName);
      double[,] rXrYrZ { get; }
      double[] rXrYrZFixedAngle { get; }
      double[,] rYrXrZ { get; }
      double[,] rZrYrX { get; }
      void StoreToCNC(string cnc_address);
      void ToFile(string FileOut, bool RobotFile);
      string[] ToString();
      void WriteAttribute(string AttribName, double Value);
      void WriteAxisAttribute(string axName, string attributeName, double attributeValue);
      bool WriteAxisPosition(string name, double MachinePos);
      double X { get; }
      double Y { get; }
      double Z { get; }
    }
  }
}
