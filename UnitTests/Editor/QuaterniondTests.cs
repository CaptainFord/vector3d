using System;
using UnityEngine;
using System.Collections.Generic;
using NUnit.Framework;
using System.Text;

namespace UnityTest {

	[TestFixture]
	[Category ("Sample Tests")]
	public class QuaterniondTests
	{
		//	This is literally 1/2^23. Multiply it by the value and that's
		//	pretty close to the precision. Rounding errors, of course, 
		//	may jitter it in any number of ways. And if the calculation
		//	involves multiple instances of rounding? The jitter may be bigger.
		const double floatPrecision = 1.19209289550781E-07;
//		const double defaultToleranceFactor = 1.19209289550781E-07;

		const int numberOfTestItems = 12;
		const int expandedNumberOfTestItems = 64;

		TestItemSet[] testItemSets;
		Dictionary<string,int[]> combos = new Dictionary<string, int[]>();

		public struct TestItemSet {
			public Quaternion fq0, fq1;
			public Quaterniond dq0, dq1;
			public Vector3 fv0, fv1;
			public Vector3d dv0, dv1;
			public float f0;
			public double d0;
		}


		public QuaterniondTests ()
		{
			System.Random rand = new System.Random("large traffic cones".GetHashCode());
			rand = new System.Random();
			testItemSets = new TestItemSet[expandedNumberOfTestItems];
			for(int i=0; i<expandedNumberOfTestItems; ++i){
				testItemSets[i] = GenerateTestItemSet(i, rand);
			}
		}

		TestItemSet GenerateTestItemSet(int index, System.Random rand){
			TestItemSet set = new TestItemSet();
			set.fq0 = GenerateRandomQuaternion(rand);
			set.fq1 = GenerateRandomQuaternion(rand);
			set.fv0 = GenerateRandomVector3(rand);
			set.fv1 = GenerateRandomVector3(rand);
			set.f0 = (float)rand.NextDouble();
			
			set.dq0 = new Quaterniond(set.fq0);
			set.dq1 = new Quaterniond(set.fq1);
			set.dv0 = new Vector3d(set.fv0);
			set.dv1 = new Vector3d(set.fv1);
			set.d0 = set.f0;

//			TestNormalizationOfSet(set);

			return set;
		}

		void TestNormalizationOfSet(TestItemSet set){
			
			if(!IsNormalized(set.fq0, floatPrecision)){
				Debug.LogWarning("Not normalized: " + set.fq0 + " = " + GetSumOfSquares(set.fq0));
			}
			if(!IsNormalized(set.fq1, floatPrecision)){
				Debug.LogWarning("Not normalized: " + set.fq1 + " = " + GetSumOfSquares(set.fq1));
			}
			if(!IsNormalized(set.dq0, floatPrecision)){
				Debug.LogWarning("Not normalized: " + set.dq0 + " = " + GetSumOfSquares(set.dq0));
			}
			if(!IsNormalized(set.dq1, floatPrecision)){
				Debug.LogWarning("Not normalized: " + set.dq1 + " = " + GetSumOfSquares(set.dq1));
			}
		}

		Quaternion GenerateRandomQuaternion(System.Random rand){
			Quaternion q = new Quaternion();
			for(int i=0; i<4; ++i){
				q[i] = (float)(2.0 * rand.NextDouble() - 1.0);	//	I want a range of negative and positive values.
			}
			Normalize (ref q);	//	Not sure if this will actually be valid or not. If it produces errors, it'll be funny at least.
			return q;
		}
		Vector3 GenerateRandomVector3(System.Random rand){
			Vector3 v3 = new Vector3();
			for(int i=0; i<3; ++i){
				//	I want a range of negative and positive values, with a hilariously large range of magnitudes.
				v3[i] = (float)(1000.0 * (0.1 + rand.NextDouble()) * (2.0 * rand.NextDouble() - 1.0));	
			}
			return v3;
		}

		void Normalize(ref Quaternion q){
			float sum = 0;
			for(int i=0; i<4; ++i){
				sum += q[i] * q[i];
			}
			float inverseMagnitude = 1f / Mathf.Sqrt(sum);
			for(int i=0; i<4; ++i){
				q[i] *= inverseMagnitude;
			}
		}

		bool IsNormalized(Quaternion q, double tolerance){
			return Math.Abs(1d - GetSumOfSquares(q)) < tolerance;
		}
		bool IsNormalized(Quaterniond q, double tolerance){
			return Math.Abs(1d - GetSumOfSquares(q)) < tolerance;
		}
		double GetSumOfSquares(Quaternion q){
			double sum = 0;
			for(int i=0; i<4; ++i){
				sum += ((double)q[i]) * q[i];
			}
			return sum;
		}
		double GetSumOfSquares(Quaterniond q){
			double sum = 0;
			for(int i=0; i<4; ++i){
				sum += q[i] * q[i];
			}
			return sum;
		}
		
		[Test(Description = "eulerAngles")]
		[Category ("eulerAngles")]
		public void TestEulerAngles (
				[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
				){
			TestItemSet set = testItemSets[testIndex];

			Vector3 fAngles;
			Vector3d dAngles;

			fAngles = set.fq0.eulerAngles;
			dAngles = set.dq1.eulerAngles;
			
			AssertSimilar(fAngles, dAngles);
		}

		[Test(Description = "SetFromToRotation")]
		[Category ("SetFromToRotation")]
		public void TestSetFromToRotation (
			[NUnit.Framework.Range (0,expandedNumberOfTestItems-1)] int testIndex
				){
			TestItemSet set = testItemSets[testIndex];

			Quaternion q = new Quaternion();
			Quaterniond qd = new Quaterniond();

			q.SetFromToRotation(set.fv0, set.fv1);
			qd.SetFromToRotation(set.dv0, set.dv1);

			if(testIndex == expandedNumberOfTestItems - 1){
				DisplayCombos();
			}

			AddOne(set, false);
			AssertSimilar(q, qd, "vectors: " + set.fv0 + ", " + set.fv1);
			AddOne(set, true);

			AssertSimilar(set.fv1.normalized, (q * set.fv0).normalized, 4d);
			AssertSimilar(set.dv1.normalized, (qd * set.dv0).normalized, 4d);
		}
		void DisplayCombos(){
			StringBuilder builder = new StringBuilder();
			foreach(string key in combos.Keys){
				int[] nums = combos[key];
				builder.Append(key + ": " + nums[0] + "/" + nums[1] + "\n");
			}
			Debug.Log (builder.ToString());

		}
		void AddOne(TestItemSet set, bool passed){
			string key = "";
			for(int i=0; i<3; ++i){
				if(set.fv0[i] < 0){
					key += "-";
				} else {
					key += "+";
				}
			}
			for(int i=0; i<3; ++i){
				if(set.fv1[i] < 0){
					key += "-";
				} else {
					key += "+";
				}
			}
			if(!combos.ContainsKey(key)){
				combos.Add (key, new int[2]);
			}
			int[] nums = combos[key];
			if(passed){
				++nums[0];
			} else {
				++nums[1];
			}
		}

		[Test(Description = "SetLookRotation")]
		[Category ("SetLookRotation")]
		public void TestSetLookRotationOneArgument (
				[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
		    	){
			TestItemSet set = testItemSets[testIndex];
			
			Quaternion q = new Quaternion();
			Quaterniond qd = new Quaterniond();
			
			q.SetLookRotation(set.fv0);
			qd.SetLookRotation(set.dv0);
			
			AssertSimilar(q, qd);
		}

		[Test(Description = "SetLookRotation")]
		[Category ("SetLookRotation")]
		public void TestSetLookRotationTwoArgument (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];
			Quaternion q = new Quaternion();
			Quaterniond qd = new Quaterniond();
			
			q.SetLookRotation(set.fv0, set.fv1);
			qd.SetLookRotation(set.dv0, set.dv1);
			
			AssertSimilar(q, qd);
		}


		[Test(Description = "ToAngleAxis")]
		[Category ("ToAngleAxis")]
		public void TestToAngleAxis (
				[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];

			Quaternion q = set.fq0;
			Quaterniond dq = set.dq0;

			float fAngle;
			Vector3 fAxis;
			
			double dAngle;
			Vector3d dAxis;
			
			q.ToAngleAxis(out fAngle, out fAxis);
			dq.ToAngleAxis(out dAngle, out dAxis);


			//	Unlike almost every other assertion, this one is based on degrees, and isn't a value between 0 and 1
			//	So that's why I need a basedOn value in the hundreds. Although, based on how close the Acos results are
			//	(they're identical to all common decimal places), I may only need to use the larger of the two angles as the basis.

//			string message = "f=" + fAngle + ";d=" + dAngle + " : " + (dAngle * Mathf.Rad2Deg / Mathd.Rad2Deg) + ";factors="
//				+ ((dAngle - fAngle) / floatPrecision) + " : " + (((dAngle * Mathf.Rad2Deg / Mathd.Rad2Deg) - fAngle) / floatPrecision + "\n");
//			message += "f.aCos=" + Mathf.Acos(set.fq0.w) + ";d.aCos=" + Mathd.Acos(set.fq0.w) + "\n";
//			message += "Angles";
			AssertSimilar(fAngle, dAngle, 1.01d * Mathd.Max(Mathd.Abs(fAngle), Mathd.Abs(dAngle)));
			AssertSimilar(fAxis, dAxis, 1.5d);
		}

		//	Static methods

		
		[Test(Description = "Static:Angle")]
		[Category ("Static:Angle")]
		public void TestStaticAngle (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];
			
			float fAngle = Quaternion.Angle(set.fq0, set.fq1);
			double dAngle = Quaterniond.Angle(set.dq0, set.dq1);
			
			AssertSimilar(fAngle, dAngle, 1.01d * Mathd.Max(Mathd.Abs(fAngle), Mathd.Abs(dAngle)));
		}

		[Test(Description = "Static:Dot")]
		[Category ("Static:Dot")]
		public void TestStaticDot (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];
			
			float fProduct = Quaternion.Dot(set.fq0, set.fq1);
			double dProduct = Quaterniond.Dot(set.dq0, set.dq1);
			
			AssertSimilar(fProduct, dProduct);
		}

		[Test(Description = "Static:Lerp")]
		[Category ("Static:Lerp")]
		public void TestStaticLerp (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];
			
			Quaternion fLerp = Quaternion.Lerp(set.fq0, set.fq1, set.f0);
			Quaterniond dLerp = Quaterniond.Lerp(set.dq0, set.dq1, set.d0);
			
			AssertSimilar(fLerp, dLerp);
		}
		
		[Test(Description = "Static:LookRotation(v3)")]
		[Category ("Static:LookRotation(v3)")]
		public void TestStaticLookRotationOneArg (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];

			Quaternion fResult = Quaternion.LookRotation(set.fv0);
			Quaterniond dResult = Quaterniond.LookRotation(set.dv0);
			
			AssertSimilar(fResult, dResult);
		}

		[Test(Description = "Static:LookRotation(v3; v3)")]
		[Category ("Static:LookRotation(v3; v3)")]
		public void TestStaticLookRotationTwoArg (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];
			
			Quaternion fResult = Quaternion.LookRotation(set.fv0, set.fv1);
			Quaterniond dResult = Quaterniond.LookRotation(set.dv0, set.dv1);
			
			AssertSimilar(fResult, dResult);
		}

		[Test(Description = "Static:RotateTowards")]
		[Category ("Static:RotateTowards")]
		public void TestStaticRotateTowards (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];
			
			Quaternion fResult = Quaternion.RotateTowards(set.fq0, set.fq1, set.f0);
			Quaterniond dResult = Quaterniond.RotateTowards(set.dq0, set.dq1, set.d0);
			
			AssertSimilar(fResult, dResult);
		}

		[Test(Description = "Static:Slerp")]
		[Category ("Static:Slerp")]
		public void TestStaticSlerp (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];
			
			Quaternion fSlerp = Quaternion.Slerp(set.fq0, set.fq1, set.f0);
			Quaterniond dSlerp = Quaterniond.Slerp(set.dq0, set.dq1, set.d0);
			
			AssertSimilar(fSlerp, dSlerp);
		}

		[Test(Description = "Identity")]
		[Category ("Static:Identity")]
		public void TestIdentity (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];
			
			AssertSimilar(Quaternion.identity, Quaterniond.identity);

			//	This should throw an exception, since I haven't implemented this operation yet. What the fuck?
			Vector3 fv = Quaternion.identity * set.fv0;
			Vector3d dv = Quaterniond.identity * set.dv0;
			
			AssertSimilar(fv, dv);
		}

		[Test]
		public void TestMultiplyQuaternions (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];

			Quaternion fq = set.fq0 * set.fq1;
			Quaterniond dq = set.dq0 * set.dq1;

			AssertSimilar(fq, dq);

			//	Testing automatic conversion. Should be identical, 
			//	since this procedure is strictly identical to the above.
			Quaterniond dq2 = set.dq0 * set.fq1;
			Assert.AreEqual(dq, dq2);
		}

		[Test]
		public void TestMultiplyVectors (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];
			
			Vector3 fv = set.fq0 * set.fv0;
			Vector3d dv = set.dq0 * set.dv0;

			AssertSimilar(fv, dv, dv.magnitude);
		}

		[Test]
		public void TestEquals (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			//	TODO Implement Quaterniond.Equals test
			//	This is actually really difficult to test.
			Assert.Fail ("Test not implemented");
		}

		[Test]
		public void TestStaticFromToRotation (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];
			
			Quaternion fResult = Quaternion.FromToRotation(set.fv0, set.fv1);
			Quaterniond dResult = Quaterniond.FromToRotation(set.dv0, set.dv1);
			
			AssertSimilar(fResult, dResult);
		}

		[Test]
		public void TestStaticAngleAxis (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];
			
			Quaternion fResult = Quaternion.AngleAxis(set.f0 * 360f, set.fv0.normalized);
			Quaterniond dResult = Quaterniond.AngleAxis(set.d0 * 360d, set.dv0.normalized);
			
			AssertSimilar(fResult, dResult, 5.0d);
		}
		[Test]
		public void TestStaticEulerVectorArg (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];
			
			Quaternion fResult = Quaternion.Euler(set.fv0);
			Quaterniond dResult = Quaterniond.Euler(set.dv0);
			
			AssertSimilar(fResult, dResult, 10.0d);
		}
		[Test]
		public void TestStaticEulerFloatArgs (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];
			
			Quaternion fResult = Quaternion.Euler(set.fv0.x, set.fv0.y, set.fv0.z);
			Quaterniond dResult = Quaterniond.Euler(set.dv0.x, set.dv0.y, set.dv0.z);
			
			AssertSimilar(fResult, dResult, 10.0d);
		}

//		[Test]	//	SPECIAL DEBUG TEST (IT WORKED!)
//		public void TestStaticEulerSpecial (
//			[NUnit.Framework.Values (false, true)] bool axis0,
//			[NUnit.Framework.Values (false, true)] bool axis1,
//			[NUnit.Framework.Values (false, true)] bool axis2,
//			[NUnit.Framework.Range (0,2)] int first,
//			[NUnit.Framework.Range (0,1)] int second
//			){
//			int[] order = {first, second >= first ? second + 1 : second, 0};
//			order[2] = 3 - order[0] - order[1];
//
//			TestItemSet set = testItemSets[0];
//			
//			Quaternion fResult = Quaternion.Euler(set.fv0);
//			Quaterniond dResult = Quaterniond.EulerSpecial(set.dv0.x, set.dv0.y, set.dv0.z,
//			                                               axis0, axis1, axis2, order);
//			
//			AssertSimilar(fResult, dResult, 100.0d);
//		}	//	END SPECIAL

		[Test]
		public void TestStaticInverse (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];
			
			Quaternion fResult = Quaternion.Inverse(set.fq0);
			Quaterniond dResult = Quaterniond.Inverse(set.dq0);
			
			AssertSimilar(fResult, dResult);
		}

		[Test]
		public void TestIndexing (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];

			for(int i=0; i<4; ++i){
				AssertSimilar(set.fq0[i], set.dq0[i], 0);
				AssertSimilar(set.fq1[i], set.dq1[i], 0);
			}
		}

		//	********************************
		//	Utility functions used for tests
		//	********************************

		void AssertSimilar(float f, double d){
			AssertSimilar(f, d, "Values");
		}

		void AssertSimilar(float f, double d, double toleranceBasedOn){
			AssertSimilar(f, d, "Values", toleranceBasedOn);
		}
		void AssertSimilar(float f, double d, String valueName){
			AssertSimilar(f, d, valueName, 1d);
		}
		void AssertSimilar(double f, double d, string valueName, double toleranceBasedOn){
			double difference = Mathd.Abs(d - f);
			double tolerance = Mathd.Abs(toleranceBasedOn * floatPrecision);
			string debugString = "f=" + f + ", \td=" + d + "\n\t"
				+ "diff=" + difference + ", \ttolerance=" + tolerance + "\n\t"
					+ "factors= " + (tolerance == 0.0d ? "&#x221E" : (difference / tolerance).ToString())
					+ " \t/\t " + (difference == 0.0d ? "&#x221E" : (tolerance / difference).ToString());
			if(difference > tolerance){
				Assert.Fail(valueName + " outside of tolerance\n\t" + debugString, f, d);
			} else {
//				if(difference != 0){
//					double passFactor = Mathd.Min ((tolerance / floatPrecision), (tolerance / difference));
//					if(passFactor >= 10.0d){
//						Debug.Log (valueName + " inside of tolerance (" + debugString + ")");
//					}
//				}
//				Assert.Pass(valueName + " inside of tolerance", f, d);
			}
		}

		void AssertSimilar(Vector3 f, Vector3 d){
			AssertSimilar(f, d, 1);
		}
		void AssertSimilar(Quaternion f, Quaternion d){
			AssertSimilar(f, d, 1);
		}
		void AssertSimilar(Vector3 f, Vector3d d){
			AssertSimilar(f, d, 1);
		}
		void AssertSimilar(Quaternion f, Quaterniond d){
			AssertSimilar(f, d, 1);
		}
		void AssertSimilar(Vector3d f, Vector3d d){
			AssertSimilar(f, d, 1);
		}
		void AssertSimilar(Quaterniond f, Quaterniond d){
			AssertSimilar(f, d, 1);
		}

		void AssertSimilar(Vector3 f, Vector3 d, double toleranceBasedOn){
			AssertSimilar(f, d, toleranceBasedOn, "");
		}
		void AssertSimilar(Quaternion f, Quaternion d, double toleranceBasedOn){
			AssertSimilar(f, d, toleranceBasedOn, "");
		}
		void AssertSimilar(Vector3 f, Vector3d d, double toleranceBasedOn){
			AssertSimilar(f, d, toleranceBasedOn, "");
		}
		void AssertSimilar(Quaternion f, Quaterniond d, double toleranceBasedOn){
			AssertSimilar(f, d, toleranceBasedOn, "");
		}
		void AssertSimilar(Vector3d f, Vector3d d, double toleranceBasedOn){
			AssertSimilar(f, d, toleranceBasedOn, "");
		}
		void AssertSimilar(Quaterniond f, Quaterniond d, double toleranceBasedOn){
			AssertSimilar(f, d, toleranceBasedOn, "");
		}

		void AssertSimilar(Vector3 f, Vector3 d, string additionalInfo){
			AssertSimilar(f, d, 1d, additionalInfo);
		}
		void AssertSimilar(Quaternion f, Quaternion d, string additionalInfo){
			AssertSimilar(f, d, 1d, additionalInfo);
		}
		void AssertSimilar(Vector3 f, Vector3d d, string additionalInfo){
			AssertSimilar(f, d, 1d, additionalInfo);
		}
		void AssertSimilar(Quaternion f, Quaterniond d, string additionalInfo){
			AssertSimilar(f, d, 1d, additionalInfo);
		}
		void AssertSimilar(Vector3d f, Vector3d d, string additionalInfo){
			AssertSimilar(f, d, 1d, additionalInfo);
		}
		void AssertSimilar(Quaterniond f, Quaterniond d, string additionalInfo){
			AssertSimilar(f, d, 1d, additionalInfo);
		}

		void AssertSimilar(Vector3 f, Vector3 d, double toleranceBasedOn, string additionalInfo){
			string inputs = (additionalInfo.Length > 1 ? additionalInfo + "\n" : "") + "f: " + f.ToString("G5") + "  d: " + d.ToString("G5") + "\n";
			AssertSimilar (f.x, d.x, inputs + "x", toleranceBasedOn);
			AssertSimilar (f.y, d.y, inputs + "y", toleranceBasedOn);
			AssertSimilar (f.z, d.z, inputs + "z", toleranceBasedOn);
		}
		void AssertSimilar(Quaternion f, Quaternion d, double toleranceBasedOn, string additionalInfo){
			string inputs = (additionalInfo.Length > 1 ? additionalInfo + "\n" : "") + "f: " + f.ToString("G5") + "  d: " + d.ToString("G5") + "\n";
			AssertSimilar (f.x, d.x, inputs + "x", toleranceBasedOn);
			AssertSimilar (f.y, d.y, inputs + "y", toleranceBasedOn);
			AssertSimilar (f.z, d.z, inputs + "z", toleranceBasedOn);
			AssertSimilar (f.w, d.w, inputs + "w", toleranceBasedOn);
		}

		void AssertSimilar(Vector3 f, Vector3d d, double toleranceBasedOn, string additionalInfo){
			string inputs = (additionalInfo.Length > 1 ? additionalInfo + "\n" : "") + "f: " + f.ToString("G5") + "  d: " + d.ToString("G5") + "\n";
			AssertSimilar (f.x, d.x, inputs + "x", toleranceBasedOn);
			AssertSimilar (f.y, d.y, inputs + "y", toleranceBasedOn);
			AssertSimilar (f.z, d.z, inputs + "z", toleranceBasedOn);
		}
		void AssertSimilar(Quaternion f, Quaterniond d, double toleranceBasedOn, string additionalInfo){
			string inputs = (additionalInfo.Length > 1 ? additionalInfo + "\n" : "") + "f: " + f.ToString("G5") + "  d: " + d.ToString("G5") + "\n";
			AssertSimilar (f.x, d.x, inputs + "x", toleranceBasedOn);
			AssertSimilar (f.y, d.y, inputs + "y", toleranceBasedOn);
			AssertSimilar (f.z, d.z, inputs + "z", toleranceBasedOn);
			AssertSimilar (f.w, d.w, inputs + "w", toleranceBasedOn);
		}

		void AssertSimilar(Vector3d f, Vector3d d, double toleranceBasedOn, string additionalInfo){
			string inputs = (additionalInfo.Length > 1 ? additionalInfo + "\n" : "") + "f: " + f.ToString("G5") + "  d: " + d.ToString("G5") + "\n";
			AssertSimilar (f.x, d.x, inputs + "x", toleranceBasedOn);
			AssertSimilar (f.y, d.y, inputs + "y", toleranceBasedOn);
			AssertSimilar (f.z, d.z, inputs + "z", toleranceBasedOn);
		}
		void AssertSimilar(Quaterniond f, Quaterniond d, double toleranceBasedOn, string additionalInfo){
			string inputs = (additionalInfo.Length > 1 ? additionalInfo + "\n" : "") + "f: " + f.ToString("G5") + "  d: " + d.ToString("G5") + "\n";
			AssertSimilar (f.x, d.x, inputs + "x", toleranceBasedOn);
			AssertSimilar (f.y, d.y, inputs + "y", toleranceBasedOn);
			AssertSimilar (f.z, d.z, inputs + "z", toleranceBasedOn);
			AssertSimilar (f.w, d.w, inputs + "w", toleranceBasedOn);
		}
	}
}

