define([
    'app',
    'text!app/templates/detail.hbs',
    'app/util/sticky',
    'app/model/detail'
], function (App, template, Sticky) {

    //var defaults = {};

    App.View.DetailButton = App.Backbone.View.extend({
        tagName: 'div',

        template: App.Handlebars.compile(
            '<button class="btn-circle" type="button">' +
                '<span class="eai eai-{{icon}}"></span>' +
            '</button>' +
            '{{#if caption}}<span>{{caption}}</span>{{/if}}'
        ),

        events: {
            "click": 'onClick'
        },

        render: function() {
            this.$el.html(this.template(this.model.attributes));

            return this;
        },

        onClick: function(e) {
            e.preventDefault();

            this.model.get('callback')(e);
        }
    });

    App.View.Detail = App.Backbone.View.extend({
        tagName: 'div',
        id: 'detail',

        sticky: null,
        template: App.Handlebars.compile(template),

        events: {
            "click": 'onClick'
        },

        clear: function () {
            if (this.model) {
                var view = this.model.get('view');
                this.off(view, 'close');
                this.off(view, 'rendered');
                this.model.destroy();
                this.model = null;
            }
            this.clearSubViews();
        },

        setModel: function (model) {
            if (!model) {
                return;
            }

            this.clear();

            this.model = model;

            var that = this,
                $bar = this.$('.bar > .wrapper > .inner').empty(),
                $document = this.$('.document > .wrapper').empty();

            App._.forEach(['left', 'right', 'center'], function(name) {
                var button = that.model.get(name + '_button');
                if (button) {
                    var view = new App.View.DetailButton({
                        model: button,
                        attributes: {
                            'class': 'detail-nav detail-nav-' + name,
                            title: button.get('title')
                        }
                    });
                    $bar.append(view.render().$el);
                    that.pushSubView(view);
                }
            });

            var view = this.model.get('view').render();
            $document.append(view.$el);
            that.pushSubView(view);

            view.postRender();

            this.listenToOnce(view, 'close', that.close);
            this.listenTo(view, 'rendered', that.applySticky);
            that.applySticky();

            this.$('.document').scrollTop(0);

            return this.open();
        },

        render: function () {
            App.$('div[data-detail-placeholder]').replaceWith(this.$el.html(this.template({})));

            return this;
        },

        applySticky: function() {
            if (this.sticky) {
                this.sticky.disable();
                this.sticky = null;
            }
            this.sticky = new Sticky({
                element: '#detail-actions',
                parent: '#detail-document',
                top: 20,
                bottom: 20
            }).enable();
        },

        onClick: function(e) {
            if (App.$(e.target).attr('id') === 'detail-document') {
                this.close();
            }
        },

        open: function() {
            App.$('body').addClass('detail-opened');

            return this;
        },

        close: function() {
            App.$('body').removeClass('detail-opened');

            this.clear();

            return this;
        }
    });

    return new App.View.Detail().render();

    //return App.View.Detail;
});
