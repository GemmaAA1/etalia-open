define(['jquery', 'app/util/utils'], function ($, Utils) {

    var values = {
        7: {icon: 'W', label: 'Week'},
        30: {icon: '1m', label: '1 month'},
        60: {icon: '2m', label: '2 months'}
    };

    var Timespan = function (options) {
        this.config = $.extend({
            debug: false,
            element: '#timespan',
            toggle: '#toggle-timespan',
            selection: '#timespan-selection'
        }, options);

        this.$element = $(this.config.element);
        this.$toggle = $(this.config.toggle);
        this.$selection = $(this.config.selection);
    };

    Timespan.getValues = function() {
        return values;
    };

    Timespan.getValueIcon = function(value) {
        if (values.hasOwnProperty(value)) {
            return values[value].icon;
        }
        throw 'Unexpected timespan value';
    };

    Timespan.getValueLabel = function(value) {
        if (values.hasOwnProperty(value)) {
            return values[value].label;
        }
        throw 'Unexpected timespan value';
    };

    Timespan.prototype.init = function() {
        var that = this,
            $body = $('body');

        this.$toggle.on('click', function(e) {
            Utils.toggleClass(that.$element, 'opened');

            e.preventDefault();
            return false;
        });
        // Close on click out
        $(window).on('click', function(e) {
            if (0 == $(e.target).closest(that.config.element).length) {
                that.$element.removeClass('opened');
            }
        });

        this.$element.on('click', '.choices a', function(e) {
            var value = $(e.target).closest('a').data('timespan');
            that.setValue(value);

            that.$toggle.trigger('click');

            $body.trigger('etalia.control.timespan.change', {
                value: value,
                label: Timespan.getValueLabel(value)
            });

            e.preventDefault();
            return false;
        });

        return this;
    };

    Timespan.prototype.setValue = function(value) {
        this.$selection
            .html(Timespan.getValueIcon(value))
            .data('value', value);

        return this;
    };

    Timespan.prototype.getValue = function() {
        return this.$selection.data('value');
    };

    return Timespan;
});
