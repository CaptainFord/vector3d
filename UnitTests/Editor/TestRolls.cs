using System;
using UnityEngine;
using System.Collections.Generic;
using NUnit.Framework;
using System.Text;



namespace UnityTest
{
	public static class Rolls {
		public static MultipliedRolls m(params NumberRoll[] rolls){
			return mult(rolls);
		}
		public static MultipliedRolls mult(params NumberRoll[] rolls){
			return new MultipliedRolls(rolls);
		}
		public static RangeRoll r(float minimum, float maximum){
			return range(minimum, maximum);
		}
		public static RangeRoll range(float minimum, float maximum){
			return new RangeRoll(minimum, maximum);
		}
		public static AddedRolls a(params NumberRoll[] rolls){
			return add(rolls);
		}
		public static AddedRolls add(params NumberRoll[] rolls){
			return new AddedRolls(rolls);
		}
		public static SelectRoll s(params float[] options){
			return select(options);
		}
		public static SelectRoll select(params float[] options){
			return new SelectRoll(options);
		}
		public static SelectRolls s(params NumberRoll[] options){
			return select(options);
		}
		public static SelectRolls select(params NumberRoll[] options){
			return new SelectRolls(options);
		}
	}
	public interface NumberRoll {
		float Roll(System.Random rand);
		NumberRoll Clone();
	}
	public class RangeRoll : NumberRoll {
		public float minimum = 0f;
		public float maximum = 1f;
		
		public RangeRoll(){
			
		}
		
		public RangeRoll (float minimum, float maximum) {
			this.minimum = minimum;
			this.maximum = maximum;
		}
		
		float NumberRoll.Roll(System.Random rand) {
			return (float)(minimum + rand.NextDouble() * (maximum - minimum));
		}

		public NumberRoll Clone(){
			return (RangeRoll)this.MemberwiseClone();
		}
	}
	public class MultipliedRolls : NumberRoll {
		public NumberRoll[] rolls;
		
		public MultipliedRolls(params NumberRoll[] rolls){
			this.rolls = rolls;
		}
		public MultipliedRolls(MultipliedRolls copyThis){
			this.rolls = UtilFunc.CopyOf(copyThis.rolls);
			int end = rolls.Length;
			for(int i=0; i<end; ++i){
				rolls[i] = rolls[i].Clone();
			}
		}
		
		float NumberRoll.Roll(System.Random rand) {
			float product = 1f;
			foreach(NumberRoll roll in rolls){
				product *= roll.Roll(rand);
			}
			return product;
		}

		public NumberRoll Clone(){
			return new MultipliedRolls(this);
		}
	}
	public class AddedRolls : NumberRoll {
		public NumberRoll[] rolls;
		
		public AddedRolls(params NumberRoll[] rolls){
			this.rolls = rolls;
		}
		public AddedRolls(AddedRolls copyThis){
			this.rolls = UtilFunc.CopyOf(copyThis.rolls);
			int end = rolls.Length;
			for(int i=0; i<end; ++i){
				rolls[i] = rolls[i].Clone();
			}
		}
		
		float NumberRoll.Roll(System.Random rand) {
			double sum = 0;
			foreach(NumberRoll roll in rolls){
				sum += roll.Roll(rand);
			}
			return (float)sum;
		}

		public NumberRoll Clone(){
			return new AddedRolls(this);
		}
	}
	public class SelectRoll : NumberRoll {
		public float[] options;
		
		public SelectRoll(params float[] options){
			this.options = options;
		}
		public SelectRoll(SelectRoll copyThis){
			this.options = UtilFunc.CopyOf(copyThis.options);
		}
		
		float NumberRoll.Roll(System.Random rand) {
			if(options.Length == 0){
				return 0;
			}
			return options[rand.Next() % options.Length];
		}

		public NumberRoll Clone(){
			return new SelectRoll(this);
		}
	}
	public class SelectRolls : NumberRoll {
		public NumberRoll[] options;
		
		public SelectRolls(params NumberRoll[] options){
			this.options = options;
		}
		
		float NumberRoll.Roll(System.Random rand) {
			if(options.Length == 0){
				return 0;
			}
			return options[rand.Next() % options.Length].Roll(rand);
		}

		public NumberRoll Clone(){
			return new SelectRolls(this);
		}
	}
	public class NumberRollProfile {
		public NumberRoll numberRoll = new RangeRoll();
		public float negativeChance = 0f;
		
		public float extremeValueChance = 0f;
		
		public float nanChance, positiveInfinityChance, negativeInfinityChance,
				maxValueChance, minValueChance, negativeMaxValueChance,
				negativeMinValueChance;
		
		public float normalizeChance;
		
		NumberRollProfile[] descendants = null;
		List<MixedRoll> mixedRolls = new List<MixedRoll>();
		float totalRelativeChance = 1f;
		float totalExtremeValueChance;
		
		public NumberRollProfile(){
			
		}
		
		public NumberRollProfile(float min, float max){
			this.numberRoll = new RangeRoll(min, max);
		}
		
		public NumberRollProfile(NumberRoll numberRoll){
			this.numberRoll = numberRoll;
		}
		public NumberRollProfile Clone(){
			NumberRollProfile clone = (NumberRollProfile)this.MemberwiseClone();
			clone.FinishBeingCloned(this);
			return clone;
		}
		internal void FinishBeingCloned(NumberRollProfile original){
			numberRoll = numberRoll.Clone();
			if(descendants != null){
				int end = descendants.Length;
				descendants = UtilFunc.CopyOf(descendants, end);
				for(int i=0; i<end; ++i){
					if(descendants[i] != null){
						descendants[i] = descendants[i].Clone ();
					}
				}
			}
			mixedRolls = new List<MixedRoll>();
			foreach(MixedRoll roll in original.mixedRolls){
				mixedRolls.Add(roll.Clone());
			}
		}
		
		public NumberRollProfile SetNegativeChance(float negativeChance){
			this.negativeChance = negativeChance;
			return this;
		}
		
		public NumberRollProfile SetExtremeValueChance(float extremeValueChance){
			this.extremeValueChance = extremeValueChance;
			return this;
		}

		float[] GetExtremeValueChances(){
			return new float[] {nanChance, positiveInfinityChance, negativeInfinityChance,
				maxValueChance, minValueChance, negativeMaxValueChance,
				negativeMinValueChance};
		}
		static readonly float[] extremeValues = {float.NaN, float.PositiveInfinity, float.NegativeInfinity,
				float.MaxValue, float.MinValue, -float.MaxValue, -float.MinValue};
		
		public NumberRollProfile GetProfileFor(int valueIndex){
			if(descendants != null && valueIndex >=0 && valueIndex < descendants.Length){
				if(descendants[valueIndex] != null){
					return descendants[valueIndex].GetProfileFor(valueIndex);
				}
			}
			return this;
		}
		public NumberRollProfile AddOverride(int valueIndex, NumberRollProfile profile){
			if(valueIndex < 0)
				throw new Exception("Value index cannot be negative (valueIndex=" + valueIndex + ")");
			if(this.descendants == null){
				this.descendants = new NumberRollProfile[valueIndex + 1];
			} else if(this.descendants.Length < valueIndex + 1){
				this.descendants = UtilFunc.CopyOf(this.descendants, valueIndex + 1);
			}
			this.descendants[valueIndex] = profile;
			return this;
		}
		public NumberRollProfile AddMixedRoll(float relativeChance, NumberRollProfile competitor){
			if(relativeChance <= 0f){
				Debug.LogError("relativeChance should be positive (relativeChance=" + relativeChance + ";profile=" + competitor + ")");
				return this;
			}
			mixedRolls.Add (new MixedRoll(relativeChance, competitor));
			totalRelativeChance += relativeChance;
			return this;
		}
		NumberRollProfile selectFromMixedRolls(System.Random rand){
			if(mixedRolls.Count == 0){
				return this;
			}
			float roll = (float)rand.NextDouble() * totalRelativeChance;
			foreach(MixedRoll mixedRoll in mixedRolls){
				if(mixedRoll.chance < roll){
					return mixedRoll.profile;
				} else {
					roll -= mixedRoll.chance;
				}
			}
			return this;
		}
		public float Roll(int valueIndex, System.Random rand){
			return GetProfileFor(valueIndex).Roll(rand);
		}
		
		public float Roll(System.Random rand){
			NumberRollProfile actualRoller = selectFromMixedRolls(rand);
			if(actualRoller != this){
				return actualRoller.Roll(rand);
			}
			CalculateTotalExtremeValueChance();
			if(totalExtremeValueChance > 0){
				if(CheckChance(extremeValueChance, rand)){
					return RollExtremeValue(rand);
				}
			}
			float value = numberRoll.Roll(rand);
			return CheckChance(negativeChance, rand) ? -value : value;
		}
		
		bool CheckChance(float chance, System.Random rand){
			if(chance <= 0f){
				return false;
			} else if(chance >= 1f){
				return true;
			} else {
				return rand.NextDouble() < chance;
			}
		}
		bool CheckSubChance(float chance, ref float roll){
			if(roll < chance){
				return true;
			} else {
				roll -= chance;
				return false;
			}
		}
		
		void CalculateTotalExtremeValueChance(){
			totalExtremeValueChance = 0f;
			foreach(float f in GetExtremeValueChances()){
				totalExtremeValueChance += Mathf.Max(0f, f);
			}
		}
		
		public float RollExtremeValue(System.Random rand) {
			float roll = (float)(rand.NextDouble() * totalExtremeValueChance);
			float[] chances = GetExtremeValueChances();
			rand.Next();
			int end = chances.Length;
			for(int i=0; i<end; ++i){
				if(CheckSubChance(chances[i],ref roll)){
					return extremeValues[i];
				}
			}
			return extremeValues[end-1];
		}
		
		class MixedRoll {
			public float chance;
			public NumberRollProfile profile;
			
			public MixedRoll (float chance, NumberRollProfile profile)
			{
				this.chance = chance;
				this.profile = profile;
			}

			public MixedRoll Clone(){
				return new MixedRoll(chance, profile.Clone());
			}
		}
	}
}