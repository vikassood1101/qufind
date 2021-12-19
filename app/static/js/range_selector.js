$(document).ready(function(){
	$('.input-left').change(function(){
		setLeftValue($(this), '');
		$(this).closest('.multi-range-slider').siblings('.range_text_box.min_range').val($(this).val());
		let current_min_Value = $(this).val().trim();
		let current_max_Value = $(this).siblings('.input-right').val().trim();
		if ( current_min_Value != '' && current_max_Value != '' ) {
			let hiddenValue = current_min_Value + '_' + current_max_Value;
			$(this).closest('.multi-range-slider').siblings('input[type=hidden]').val(hiddenValue);
		} 
	});
	$('.input-right').change(function(){
		setRightValue($(this), '');
		$(this).closest('.multi-range-slider').siblings('.range_text_box.max_range').val($(this).val());
		let current_min_Value = $(this).siblings('.input-left').val().trim();
		let current_max_Value = $(this).val().trim();
		if ( current_min_Value != '' && current_max_Value != '' ) {
			let hiddenValue = current_min_Value + '_' + current_max_Value;
			$(this).closest('.multi-range-slider').siblings('input[type=hidden]').val(hiddenValue);
		}
	});
	
	$('.range_text_box').on('keyup',function(e) {
		if((e.which || e.keyCode) == 13) {
			range_text_enter($(this));
		}
	});

	$('.range_text_box').blur(function(){
		range_text_enter($(this));
	});
	
	$('.input-left').mouseover(function() {
		$(this).siblings('.slider').children('.thumb.left').addClass("hover");
	});
	$('.input-left').mouseout(function() {
		$(this).siblings('.slider').children('.thumb.left').removeClass("hover");
	});
	$('.input-left').mousedown(function() {
		$(this).siblings('.slider').children('.thumb.left').addClass("active");
	});
	$('.input-left').mouseup(function() {
		$(this).siblings('.slider').children('.thumb.left').removeClass("active");
	});
	
	$('.input-right').mouseover(function() {
		$(this).siblings('.slider').children('.thumb.right').addClass("hover");
	});
	$('.input-right').mouseout(function() {
		$(this).siblings('.slider').children('.thumb.right').removeClass("hover");
	});
	$('.input-right').mousedown(function() {
		$(this).siblings('.slider').children('.thumb.right').addClass("active");
	});
	$('.input-right').mouseup(function() {
		$(this).siblings('.slider').children('.thumb.right').removeClass("active");
	});
});

function range_text_enter(thisElement){
	let min_fixed_value = parseInt($(thisElement).attr('min').trim()); 
	let max_fixed_value = parseInt($(thisElement).attr('max').trim()); 
	let current_value = parseInt($(thisElement).val().trim());
	if ($(thisElement).hasClass('min_range')){
		if ( min_fixed_value <= current_value && current_value < max_fixed_value && 
			current_value < parseInt($(thisElement).siblings('.max_range').val().trim()) ) {
				$(thisElement).siblings('.multi-range-slider').children('.input-left').val(current_value);
				$(thisElement).siblings('input[type=hidden]').val($(thisElement).val() + '_' + $(thisElement).siblings('.max_range').val());
		} else {
			$(thisElement).siblings('input[type=hidden]').val(min_fixed_value.toString() + '_' + $(thisElement).siblings('.max_range').val());
			$(thisElement).val(min_fixed_value);
			$(thisElement).siblings('.multi-range-slider').children('.input-left').val(min_fixed_value);
			current_value = min_fixed_value;
		}
		setLeftValue($(thisElement).siblings('.multi-range-slider').children('.input-left'), current_value);
	} else {
		if ( min_fixed_value < current_value && current_value <= max_fixed_value &&
		current_value > parseInt($(thisElement).siblings('.min_range').val().trim()) ) {
			$(thisElement).siblings('.multi-range-slider').children('.input-right').val(current_value);
			$(thisElement).siblings('input[type=hidden]').val($(thisElement).siblings('.min_range').val() + '_' + $(thisElement).val());
		} else {
			$(thisElement).siblings('input[type=hidden]').val($(thisElement).siblings('.min_range').val() + '_' + max_fixed_value.toString());
			$(thisElement).val(max_fixed_value);
			$(thisElement).siblings('.multi-range-slider').children('.input-right').val(max_fixed_value);
			current_value = max_fixed_value;
		}
		setRightValue($(thisElement).siblings('.multi-range-slider').children('.input-right'), current_value);
	}
}

function setLeftValue(thisElement, labelValue) {
	let min = parseInt($(thisElement).attr('min'));
	let	max = parseInt($(thisElement).attr('max'));

	$(thisElement).val(Math.min(parseInt($(thisElement).val()), parseInt($(thisElement).siblings('.input-right').val()) - 1));
	let percent = ((parseInt($(thisElement).val()) - min) / (max - min)) * 100;
	
	if (labelValue == ''){
		$(thisElement).siblings('.slider').find('.thumb.left>.thumb_label').text($(thisElement).val());
	} else {
		$(thisElement).siblings('.slider').find('.thumb.left>.thumb_label').text(labelValue);
	}

	$(thisElement).siblings('.slider').children('.thumb.left').css({ left : percent + "%" });
	$(thisElement).siblings('.slider').children('.range').css({ left : percent + "%" });
}

function setRightValue(thisElement, labelValue) {
	let min = parseInt($(thisElement).attr('min'));
	let	max = parseInt($(thisElement).attr('max'));

	$(thisElement).val(Math.max(parseInt($(thisElement).val()), parseInt($(thisElement).siblings('.input-left').val()) + 1));

	if (labelValue == ''){
		$(thisElement).siblings('.slider').find('.thumb.right>.thumb_label').text($(thisElement).val());
	} else {
		$(thisElement).siblings('.slider').find('.thumb.right>.thumb_label').text(labelValue);
	}
	
	let percent = ((parseInt($(thisElement).val())  - min) / (max - min)) * 100;
    $(thisElement).siblings('.slider').children('.thumb.right').css({ right : (100 - percent) + "%" });
	$(thisElement).siblings('.slider').children('.range').css({ right : (100 - percent) + "%" });
}

