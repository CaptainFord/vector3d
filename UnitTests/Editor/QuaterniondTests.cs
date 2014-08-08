using System;
using UnityEngine;
using System.Collections.Generic;
using NUnit.Framework;
using System.Text;

namespace UnityTest {


	[TestFixture]
	[Category ("Quaternion")]
	public class QuaterniondTests : TestsCommon
	{
		//	This is literally 1/2^23. Multiply it by the value and that's
		//	pretty close to the precision. Rounding errors, of course, 
		//	may jitter it in any number of ways. And if the calculation
		//	involves multiple instances of rounding? The jitter may be bigger.

//		const double defaultToleranceFactor = 1.19209289550781E-07;

		const int numberOfTestItems = 1;
		const int expandedNumberOfTestItems = 12;
		public TestItemSetCollection testItemSets;

		IDictionary<string,ComboInfo> combos = new SortedList<string, ComboInfo>();
//		List<int> list = new List<int>();

		public struct TestItemSet {
			public Quaternion fq0, fq1;
			public Quaterniond dq0, dq1;
			public Vector3 fv0, fv1, fEuler;
			public Vector3d dv0, dv1, dEuler;
			public float f0;
			public double d0;

			public override string ToString(){
				return ToString("G5");
			}
			public string ToString(string format){
				return "q0=" + dq0.ToString(format) + " q1=" + dq1.ToString(format) 
					+ " v0=" + dv0.ToString(format) + " v1=" + dv1.ToString(format) 
					+ " d0=" + d0.ToString(format);
			}
		}

		public class TestItemSetCollection {
			TestItemSet[] testItemSets;
			QuaterniondTests parent;

			public TestItemSetCollection (QuaterniondTests parent)
			{
				this.parent = parent;
			}

			public TestItemSet this[int index] {
				get {
					if(testItemSets == null){
						testItemSets = parent.GenerateTestItemSets();
					}
					return testItemSets[index];
				}
			}
		}

		public QuaterniondTests() {
			testItemSets = new TestItemSetCollection(this);
		}

		internal TestItemSet[] GenerateTestItemSets ()
		{
			System.Random rand = new System.Random("large traffic cones".GetHashCode());
//			rand = new System.Random();
			int length = Math.Max (expandedNumberOfTestItems, numberOfTestItems);
			TestItemSet[] testItemSets = new TestItemSet[length];
			for(int i=0; i<length; ++i){
				testItemSets[i] = GenerateTestItemSet(i, rand);
			}
			return testItemSets;
		}

		TestItemSet GenerateTestItemSet(int index, System.Random rand){
			TestItemSet set = new TestItemSet();
			set.fq0 = GenerateRandomQuaternion(rand);
			set.fq1 = GenerateRandomQuaternion(rand);
			set.fv0 = GenerateRandomVector3(rand);
			set.fv1 = GenerateRandomVector3(rand);
			set.fEuler = GenerateRandomEulerAngles(rand);
			set.f0 = (float)rand.NextDouble();
			
			set.dq0 = new Quaterniond(set.fq0);
			set.dq1 = new Quaterniond(set.fq1);
			set.dv0 = new Vector3d(set.fv0);
			set.dv1 = new Vector3d(set.fv1);
			set.dEuler = new Vector3d(set.fEuler);
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

		Vector3 GenerateRandomEulerAngles(System.Random rand){
			Vector3 v3 = new Vector3();
//			v3.x = RollBoundedDouble(rand, -90f, 90f, 0.2f);
//			v3.x = RollBoundedDouble(rand, -180f, 180f, 0.2f);
//			v3.x = RollBoundedDouble(rand, -180f, 180f, 0.2f);
			v3.x = (float)(rand.NextDouble() * 180) - 90f;
			if(v3.x < 0f){
				v3.x += 360f;
			}
			v3.y = (float)(rand.NextDouble() * 360);
			v3.z = (float)(rand.NextDouble() * 360);
//			v3.y = RollBoundedDouble(rand, 0, 360f, 0.2f);
//			v3.z = RollBoundedDouble(rand, 0, 360f, 0.2f);
			return v3;

		}

		float RollBoundedDouble(System.Random rand, float min, float max, float boundaryChance){
			if(rand.NextDouble() < boundaryChance){
				switch(rand.Next() % 3){
				case 0:
					return min;
				case 1:
					return max;
				default:
					return (min + max) * 0.5f;
				}
			} else {
				return (float)(rand.NextDouble() * (max - min) + min);
			}
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
			dAngles = set.dq0.eulerAngles;
			
			AssertSimilar(fAngles, dAngles, 360d);
		}

		public struct EulerMatch : IComparable {
			public bool flippedMatrix;
			public int[] tuple;
			public string valueReorder;
			public int[] returnOrder;

			public Vector3d eulerAngles;
			public Quaterniond reconstruction;

			public Vector3d resultA;
			public Vector3d resultB;

			public double totalDiffSquared;

			public EulerMatch(int tupleNum, bool flipMatrix, int reorderNum){
				tuple = new int[4];
				tuple[0] = (tupleNum % 3) + 1;
				tuple[1] = (tupleNum / 3) % 2;
				tuple[2] = (tupleNum / 6) % 2;
				tuple[3] = (tupleNum / 12) % 2;
				this.flippedMatrix = flipMatrix;
				returnOrder = new int[3];
				returnOrder[0] = reorderNum % 3;
				returnOrder[1] = (reorderNum / 3) % 2;
				if(returnOrder[1] >= returnOrder[0]){
					++returnOrder[1];
				}
				returnOrder[2] = 3 - returnOrder[0] - returnOrder[1];
				string[] labels = {"x","y","z"};
				valueReorder = labels[returnOrder[0]]+labels[returnOrder[1]]+labels[returnOrder[2]];

				eulerAngles = new Vector3d();
				reconstruction = new Quaterniond();
				resultA = new Vector3d();
				resultB = new Vector3d();
				totalDiffSquared = 0d;
			}

			public override string ToString(){
				StringBuilder msg = new StringBuilder();
				for(int i=0; i<4; ++i){
					msg.Append(tuple[i]);
				}
				if(flippedMatrix){
					msg.Append("f");
				}
				msg.Append("-").Append(valueReorder);
				msg.Append(" = ").Append(resultA).Append(", ").Append(resultB);
				msg.Append(" => ").Append(totalDiffSquared);
				msg.Append(" (euler=" + eulerAngles.ToString("G5") + ";quat=" + reconstruction + ")");
				return msg.ToString();
			}

			public void EvaluateAngles(Quaterniond original){
				this.eulerAngles = original.EulerAngles(this.flippedMatrix, this.tuple);
				this.eulerAngles = new Vector3d(this.eulerAngles[returnOrder[0]],this.eulerAngles[returnOrder[1]],this.eulerAngles[returnOrder[2]]);
				this.reconstruction = Quaterniond.Euler(this.eulerAngles);
			}

			public void EvaluateResults(Vector3d inputA, Vector3d inputB){
				this.resultA = this.reconstruction * inputA;
				this.resultB = this.reconstruction * inputB;
			}

			public double CompareResults(Vector3d actualA, Vector3d actualB){
				this.totalDiffSquared = (actualA - this.resultA).sqrMagnitude + 
					(actualB - this.resultB).sqrMagnitude;
				return totalDiffSquared;
			}
			public int CompareTo (object obj){
				if(!(obj is EulerMatch)){
					return -1;
				}
				EulerMatch other = (EulerMatch)obj;
				return System.Math.Sign(this.totalDiffSquared - other.totalDiffSquared);
			}
		}

//		[Test(Description = "eulerAngles")]
		public void FindBestEulerMatch (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];

			Vector3d actualA = set.dq0 * set.dv0;
			Vector3d actualB = set.dq0 * set.dv1;

			EulerMatch best = new EulerMatch();

			List<EulerMatch> allAttempts = new List<EulerMatch>();

			for(int i=0; i<24*2*6; ++i){
				EulerMatch match = new EulerMatch(i, (i / 24) % 2 == 0, i / 48);
				match.EvaluateAngles(set.dq0);
				match.EvaluateResults(set.dv0, set.dv1);
				match.CompareResults(actualA, actualB);
				allAttempts.Add(match);
				if(i == 0){
					best = match;
				} else {
					if(match.totalDiffSquared < best.totalDiffSquared){
						best = match;
					}
				}
			}

			allAttempts.Sort();

			StringBuilder message = new StringBuilder();
			message.Append("Actual Results: " + actualA + ", " + actualB);
			message.Append("\nBest Match: " + best);
			message.Append("\nInput Quaternion: " + set.dq0 + 
			               "\nUnity's Euler Result: " + set.fq0.eulerAngles.ToString("G5"));
			foreach(EulerMatch match in allAttempts){
				message.Append("\n").Append(match.ToString());
			}

			Debug.Log (message.ToString());
			AssertSimilar(ModAnglesToMatch((Vector3d)set.fq0.eulerAngles, best.eulerAngles), best.eulerAngles, 2 * 360d);
		}

		[Test(Description = "eulerAngles")]
		[Category ("eulerAngles")]
		public void TestEulerAngles2 (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];

			Quaternion q = Quaternion.Euler(set.fEuler);
			Quaterniond dq = (Quaterniond)q;

			Vector3 fAngles;
			Vector3d dAngles;
			
			fAngles = q.eulerAngles;
			dAngles = dq.eulerAngles;
			
			AssertSimilar(fAngles, dAngles, 360d);
		}

		[Test]
		[Category ("eulerAngles")]
		public void TestEulerAnglesConsistency (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];
			
			Vector3 fAngles;
			Vector3d dAngles;
			
			fAngles = set.fq0.eulerAngles;
			dAngles = set.dq1.eulerAngles;

			Quaternion fresult = Quaternion.Euler(fAngles);
			Quaterniond dresult = Quaterniond.Euler(dAngles);

			AssertSimilar(fresult, set.fq0, 2d, "Unity's Quaternion");
			AssertSimilar(dresult, set.dq0);
		}

		[Test]
		[Category ("eulerAngles")]
		public void TestEulerAnglesConsistency2 (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];

			Quaternion fresult = Quaternion.Euler(set.fEuler);
			Quaterniond dresult = Quaterniond.Euler(set.dEuler);

			Vector3 fAngles = fresult.eulerAngles;
			Vector3d dAngles = dresult.eulerAngles;

			Vector3 fExpected = set.fEuler;
			Vector3d dExpected = set.dEuler;

			for(int i=0; i<3; ++i){
				fExpected[i] = ModAngleToMatch(fExpected[i], fAngles[i]);
				dExpected[i] = ModAngleToMatch(dExpected[i], dAngles[i]);
			}
			
//			AssertSimilar(fExpected, fAngles, 4*360.01d, "Unity's Quaternion");
			AssertSimilar(dExpected, dAngles, 1.5*360.01d);
		}
		float ModAngleToMatch(float value, float matchMe){
			while(value > matchMe + 180f){
				value -= 360f;
			}
			while(value < matchMe - 180f){
				value += 360f;
			}
			return value;
		}
		double ModAngleToMatch(double value, double matchMe){
			while(value > matchMe + 180d){
				value -= 360d;
			}
			while(value < matchMe - 180d){
				value += 360d;
			}
			return value;
		}

		Vector3d ModAnglesToMatch(Vector3d values, Vector3d matchMe){
			for(int i=0; i<3; ++i){
				values[i] = ModAngleToMatch(values[i], matchMe[i]);
			}
			return values;
		}

		[Test(Description = "SetFromToRotation")]
		[Category ("SetFromToRotation")]
		public void TestSetFromToRotation (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
				){
			TestItemSet set = testItemSets[testIndex];

			Quaternion q = new Quaternion();
			Quaterniond qd = new Quaterniond();

			q.SetFromToRotation(set.fv0, set.fv1);
			qd.SetFromToRotation(set.dv0, set.dv1);

//			Quaterniond qdneg = Negative(qd);
//			if(GetSumOfSquares(Subtract(q, qd)) < GetSumOfSquares(Subtract(q, qdneg))){
//				AssertSimilar(q, qd, 2d);
//			} else {
//				AssertSimilar(q, qdneg, 2d);
			//			}
			AssertSimilar(q, qd, 2d);

			AssertSimilar(set.fv1.normalized, (q * set.fv0).normalized, 4d);
			AssertSimilar(set.dv1.normalized, (qd * set.dv0).normalized, 4d);
		}




		[Test(Description = "SetLookRotation")]
		[Category ("SetLookRotation")]
		public void TestSetLookRotationOneArgument (
				[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
		    	){
			TestItemSet set = testItemSets[testIndex];
			
			Quaternion q = new Quaternion();
			Quaterniond qd = new Quaterniond();
//
//			Quaternion q2 = new Quaternion();
//			Quaterniond qd2 = new Quaterniond();

//			q2.SetFromToRotation(Vector3.forward, set.fv0);
//			qd2.SetFromToRotation(Vector3d.forward, set.dv0);
			
			q.SetLookRotation(set.fv0);
			qd.SetLookRotation(set.dv0);

//			Debug.Log ("q=" + q.ToString("G5") + " q2=" + q2.ToString("G5") + " qd=" + qd.ToString("G5") + " qd2=" + qd2.ToString("G5") +
//			           "\nq=" + Quaternion.Angle(Quaternion.identity, q) + " q2=" + Quaternion.Angle(Quaternion.identity, q2) 
//			           + " qd=" + Quaterniond.Angle(Quaterniond.identity, qd) + " qd2=" + Quaterniond.Angle(Quaterniond.identity, qd2));
//			Vector3[] vectors = new Vector3[8];
//			for(int i=0; i<8; ++i){
//				vectors[i] = new Vector3(
//					i % 2 == 0 ? -1f : 1f,
//					(i / 2) % 2 == 0 ? -1f : 1f,
//					(i / 4) == 0 ? -1f : 1f
//					);
//			}
//			StringBuilder msg = new StringBuilder();
//			foreach(Vector3 coords in vectors){
//				msg.Append(coords).Append(" -> ").Append((q * coords).ToString("G5")).Append(" => ").Append((qd * new Vector3d(coords)).ToString("G5")).AppendLine();
//			}
//			Debug.Log (msg);

			AssertSimilar (set.fv0.normalized, q * Vector3.forward, 1.5d);
			AssertSimilar (set.dv0.normalized, qd * Vector3d.forward, 1.5d);
			
			AssertSimilar(q, qd, "vector: " + set.fv0);


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

//						Vector3[] vectors = new Vector3[8];
//						for(int i=0; i<8; ++i){
//							vectors[i] = new Vector3(
//								i % 2 == 0 ? -1f : 1f,
//								(i / 2) % 2 == 0 ? -1f : 1f,
//								(i / 4) == 0 ? -1f : 1f
//								);
//						}
//						StringBuilder msg = new StringBuilder();
//						foreach(Vector3 coords in vectors){
//							msg.Append(coords).Append(" -> ").Append((q * coords).ToString("G5")).Append(" => ").Append((qd * new Vector3d(coords)).ToString("G5")).AppendLine();
//						}
//						Debug.Log (msg);

			AssertSimilar (set.fv0.normalized, q * Vector3.forward, 2.5d);
			AssertSimilar (set.dv0.normalized, qd * Vector3d.forward, 2.5d);
			
			AssertSimilar(q, qd, "vectors: " + set.fv0 + ", " + set.fv1);
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
			
			AssertSimilar(fLerp, dLerp, "q0=" + set.dq0 + " q1=" + set.dq1 + " t=" + set.f0);
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
			
			RotateTowardsTestPair(set, 360f);
			RotateTowardsTestPair(set, 90f);
			RotateTowardsTestPair(set, 30f);
			RotateTowardsTestPair(set, 15f);
			RotateTowardsTestPair(set, 4f);
			RotateTowardsTestPair(set, 0);

		}

		void RotateTowardsTestPair(TestItemSet set, float maxAngleMult){
			//	Once with the raw angle, and once with it multiplied by the input floats (might as well do both)
			RotateTowardsTest(set, maxAngleMult);
			RotateTowardsTest(set, maxAngleMult * set.f0);
		}
		void RotateTowardsTest(TestItemSet set, float angle){
			Quaternion fResult = Quaternion.RotateTowards(set.fq0, set.fq1, angle);
			Quaterniond dResult = Quaterniond.RotateTowards(set.dq0, set.dq1, angle);
			AssertSimilar(fResult, dResult, 1.5d);
		}

		[Test(Description = "Static:Slerp")]
		[Category ("Static:Slerp")]
		public void TestStaticSlerp (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];
			
			Quaternion fSlerp = Quaternion.Slerp(set.fq0, set.fq1, set.f0);
			Quaterniond dSlerp = Quaterniond.Slerp(set.dq0, set.dq1, set.d0);
			
			AssertSimilar(fSlerp, dSlerp, "q0=" + set.dq0 + " q1=" + set.dq1 + " t=" + set.f0);
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
			TestItemSet set = testItemSets[testIndex];
			string msg = set.ToString();

			EqualsTest(set.fq0, set.fq0, set.dq0, set.dq0, true, msg);
			EqualsTest(set.fq0, set.fq1, set.dq0, set.dq1, false, msg);
		}

		void EqualsTest(Quaternion q0, Quaternion q1, Quaterniond dq0, Quaterniond dq1,
		                bool sameInstance,
		                string msg) {
			if(sameInstance){
				if(!(q0.Equals(q1) && q0 == q1)){
					Assert.Fail("Unity's Quaternion fails on self-equivalence!\n" + msg);
				}
				Assert.True(dq0 == dq1, msg);
				Assert.True(dq0.Equals(dq1), msg);
			} else {
				if((q0.Equals(q1)) != (q0 == q1)){
					Assert.Fail("Unity's Quaternion result of == does not match .Equals()!\n" + msg);
				}
				Assert.AreEqual(q0.Equals(q1),dq0.Equals(dq1), msg);
				Assert.AreEqual(q0 == q1,dq0 == dq1, msg);
			}
		}

		[Test]
		public void TestApproximately (
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

		protected override void AssertSimilar(Quaternion f, Quaterniond d, double toleranceBasedOn, string additionalInfo){
//			AddCombo(f, d, false);
			AssertSimilar(new Quaterniond(f), d, toleranceBasedOn, additionalInfo);
//			AddCombo(f, d, true);
		}


		int numComboEntries = 0;
		int numMatchingPattern = 0;
		List<string> patternExceptions = new List<string>();

		void AddCombo(Quaternion f, Quaterniond d, bool passed){
			string mainKey = ComboKeyOf(d);
			string otherKey = ComboKeyOf(f);
			mainKey += GetNameOfHighestMagnitudeCoordinate(d);
			otherKey += GetNameOfHighestMagnitudeCoordinate(d);
			ComboInfo info = GetOrCreateComboInfo(mainKey);
			if(passed){
				++info.passed;
			} else {
				++info.count;
				if(mainKey.Equals(otherKey)){
					++info.matchesSignsOfDouble;
//					info.matchingValues.Add(d);
				} else {
//					info.mismatchingValues.Add(d);
				}

				++numComboEntries;
				if(MatchesTestPattern(f)){
					++numMatchingPattern;
				} else {
//					patternExceptions.Add (f.ToString("G5"));
					if(MatchesSecondaryTestPattern(f)){
						patternExceptions.Add ((f.x / f.w).ToString());
					} else {
						patternExceptions.Add (f.ToString("G5"));
					}
				}
			}
		}

		string GetNameOfHighestMagnitudeCoordinate(Quaterniond q){
			int maxIndex = 0;
			double maxValue = Mathd.Abs(q[0]);
			for(int i=1; i<4; ++i){
				double value = Mathd.Abs(q[i]);
				if(value > maxValue){
					maxIndex = i;
					maxValue = value;
				}
			}
			return "xyzw".Substring(maxIndex, 1);
		}

		bool MatchesSecondaryTestPattern(Quaternion f){
			if(f.x >= 0f || f.w < 0f){
				return false;
			}
//			float x = Mathf.Abs(f.x);
//			float y = Mathf.Abs(f.y);
//			float z = Mathf.Abs(f.z);
//			float w = Mathf.Abs(f.w);
			return false;
		}

		bool MatchesTestPattern(Quaternion f){
			//	My suspicion is that unity's quaternions always have the element with the highest magnitude as positive
			bool highestIsPositive = f[0] >= 0f;
			bool highestIsAmbiguous = false;
			float maxValue = Mathf.Abs(f[0]);
			for(int i=1; i<4; ++i){
				float value = Mathf.Abs(f[i]);
				bool isPositive = value >= 0f;
				if(value > maxValue){
					maxValue = value;
					highestIsPositive = isPositive;
					highestIsAmbiguous = false;
				} else if(value == maxValue){
					if(isPositive != highestIsPositive){
						highestIsAmbiguous = true;
					}
				}
			}
			if(highestIsAmbiguous){
				Debug.LogWarning("Tie for max value's sign: " + f.ToString("G9"));
				return true;
			}
			return highestIsPositive;
		}

		public class ComboInfo {
			public string key;
			public int passed;
			public int count;
			public int matchesSignsOfDouble;

			public List<Quaternion> matchingValues = new List<Quaternion>();
			public List<Quaternion> mismatchingValues = new List<Quaternion>();

			public override string ToString ()
			{
				return key + "(" + passed + "," + matchesSignsOfDouble + "/" + count + ")";
			}
		}

		ComboInfo GetOrCreateComboInfo(string key){
			if(!combos.ContainsKey(key)){
				ComboInfo info = new ComboInfo();
				info.key = key;
				combos.Add (key, info);
			}
			return combos[key];
		}

		string ComboKeyOf(Quaternion q){
			StringBuilder b = new StringBuilder();
			for(int i=0; i<4; ++i){
				b.Append (q[i] < 0f ? "-" : (q[i] == 0f ? "0" : "+"));
			}
			return b.ToString();
		}
		string ComboKeyOf(Quaterniond q){
			StringBuilder b = new StringBuilder();
			for(int i=0; i<4; ++i){
				b.Append (q[i] < 0d ? "-" : (q[i] == 0d ? "0" : "+"));
			}
			return b.ToString();
		}

		[TestFixtureTearDown]
		public void LogComboInfo(){
			if(combos.Count > 0){
				string message = BuildComboInfoMessage();
				Debug.Log (message);
			}
		}
//		[TearDown]
//		public static void LogComboInfo(QuaterniondTests tests){
//			string message = tests.BuildComboInfoMessage();
//			Debug.Log (message);
//		}

		string BuildComboInfoMessage(){
			StringBuilder msg = new StringBuilder();
			IList<ComboInfo> allPassed = new List<ComboInfo>();
			IList<ComboInfo> somePassed = new List<ComboInfo>();
			IList<ComboInfo> nonePassed = new List<ComboInfo>();

//			msg.Append("Pattern Matches: " + numMatchingPattern + "/" + numComboEntries);
//			AppendItems(msg, patternExceptions);
//			msg.AppendLine();

			foreach(string key in combos.Keys){
				ComboInfo info = combos[key];
				if(info.passed == info.count){
					allPassed.Add (info);
				} else if(info.passed == 0){
					nonePassed.Add(info);
				} else {
					if(info.passed < 0 || info.passed > info.count){
						Debug.LogError("Combo Info has a bizarre number of passes: " + info);
					}
					somePassed.Add(info);
				}
			}
			msg.Append("100% Passed[" + allPassed.Count + "]: (");
			appendItemNames(msg, allPassed);
			msg.Append(") ");

			msg.Append("100% Failed[" + nonePassed.Count + "]: (");
			appendItemNames(msg, nonePassed);
			msg.Append(") ");

			msg.Append("Others[" + somePassed.Count + "]: (");
			appendItemNames(msg, somePassed);
			msg.Append(")");
			foreach(ComboInfo info in combos.Values){
				msg.AppendLine();
				msg.Append(info.key).Append(" ").Append(info.passed).Append("/").Append(info.count)
						.Append(" ").Append(info.matchesSignsOfDouble).Append("/").Append(info.count);
//				msg.Append("  ");
//				AppendInfoValues(msg, info);
			}

			return msg.ToString();
		}

		void AppendInfoValues(StringBuilder msg, ComboInfo info){
			msg.AppendLine().Append("\tMatches[").Append(info.matchingValues.Count).Append("]: ");
			if(info.matchingValues.Count > 0){
				foreach(Quaternion q in info.matchingValues){
					msg.Append(q.ToString("G3")).Append(", ");
				}
				msg.Remove(msg.Length - 2, 2);
			}
			msg.AppendLine().Append("\tMismatches[").Append(info.mismatchingValues.Count).Append("]: ");
			if(info.matchingValues.Count > 0){
				foreach(Quaternion q in info.mismatchingValues){
					msg.Append(q.ToString("G3")).Append(", ");
				}
				msg.Remove(msg.Length - 2, 2);
			}
		}

		void appendItemNames(StringBuilder builder, ICollection<ComboInfo> items){
			if(items.Count < 1){
				return;
			}
			const string separator = ", ";
			foreach(ComboInfo item in items){
				builder.Append(item.key).Append(separator);
			}
			builder.Remove(builder.Length - separator.Length, separator.Length);
		}

		void AppendItems(StringBuilder builder, ICollection<Quaternion> items){
			if(items.Count < 1){
				return;
			}
			const string separator = ", ";
			foreach(Quaternion item in items){
				builder.Append(item.ToString("G5")).Append(separator);
			}
			builder.Remove(builder.Length - separator.Length, separator.Length);
		}

		void AppendItems<T>(StringBuilder builder, ICollection<T> items){
			if(items.Count < 1){
				return;
			}
			const string separator = ", ";
			foreach(T item in items){
				builder.Append(item).Append(separator);
			}
			builder.Remove(builder.Length - separator.Length, separator.Length);
		}

		/*
		 * Results from TestSetFromToRotation with 1024 tests using the seed "large traffic cones".GetHashCode()
		 * 100% Passed[8]: (---+, +--+, -+-+, --++, ++-+, +-++, -+++, ++++) 100% Failed[7]: (--+-, -+--, +---, +-+-, -++-, ++--, +++-) Others[0]: ()
			---+ 87/87 87/87
			--+- 0/11 0/11
			-+-- 0/8 0/8
			+--- 0/8 0/8
			+--+ 105/105 105/105
			-+-+ 119/119 119/119
			--++ 109/109 109/109
			+-+- 0/18 0/18
			-++- 0/21 0/21
			++-- 0/19 0/19
			++-+ 113/113 113/113
			+-++ 117/117 117/117
			-+++ 128/128 128/128
			+++- 0/36 0/36
			++++ 125/125 125/125

			//	I was RIGHT! It always keeps the real portion positive. It makes comparison so much easier.
			//	That needs to be in the Normalize method, and I also need a new method called NormalizeSign()

			//	Oh wait, no, it's the opposite. Mine always had a positive real value. Unity's version occasionally doesn't.
			//	...crap. Why? Is there anything meaningful about it? Or is it the result of some optimization I'm unaware of?

				In attempting to identify a pattern, I've concluded that there really isn't anything dependable.
				I just need to know the algorithm. And I don't think it ultimately matters.
		*/
	}
}

